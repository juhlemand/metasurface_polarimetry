/****************************************************************************

Thorlabs TXP PAX57xx single measurement demo application

This demo program connects to an instrument, configures it and takes
single measurment results in a loop. Finally it closes the connection.

Source file 'SingleMeas.c'

Date:          Jul-22-2005
Software-Nr:   N/A
Version:       1.1
Copyright:     Copyright(c) 2005, Thorlabs GmbH (www.thorlabs.com)
Author:        Michael Biebl (mbiebl@thorlabs.com)

Changelog:     Mar-01-2005 -> V1.0
Jul-22-2005 -> adjustable wavelength -> V1.1


Disclaimer:

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "TXP_DRV_PAX.h"
//#include "..\..\..\include\TXP_DRV_PAX.h"


/*---------------------------------------------------------------------------
Macros
---------------------------------------------------------------------------*/
#define DEF_IP_ADDR           "localhost"    // instrument IP address or hostname ('localhost' for own PC)
#define DEF_IP_PORT           2402           // TXP standard IP port number
#define DEF_SLOT              0              // zero based slot number
#define DEF_COMM_TIMEOUT      30000          // communication timeout
#define DEF_BOOT_TIMEOUT      30000          // wait for card to finish booting 

#define DEF_WAVELENGTH        532      		// operating wavelength 1550.0nm

#define BASIC_MEAS_SPEED      200.0          // measuremnet speed (2* motor speed) per second
#define PERIODS_PER_BUFFER    2              // take 2 half revolutions per measuremnt buffer
#define BUFFERS_PER_FFT       10             // take 10 buffers per transformation
#define FFT_PER_RESULT        2              // average 2 transformations per result

#define NUM_MEASUREMENTS      10             // number of single measurements to take


/*---------------------------------------------------------------------------
Prototypes
---------------------------------------------------------------------------*/
void errorMessage(ViSession, ViStatus);
void statusMessage(ViUInt32);
void waitKeypress(void);
void delay1sec(void);
void usage(char*);



/*---------------------------------------------------------------------------
Main
---------------------------------------------------------------------------*/
int main(int argc, char** argv)
{
	char           ip_addr[60] = DEF_IP_ADDR;
	unsigned short ip_port = DEF_IP_PORT;
	unsigned int   slot = DEF_SLOT;
	unsigned int   comm_tmo = DEF_COMM_TIMEOUT;
	unsigned int   boot_tmo = DEF_BOOT_TIMEOUT;
	unsigned int   wavelength = DEF_WAVELENGTH;

	char           resName[80];   // buffer to build resource name
	ViSession      instr;         // instruenmt handle
	ViUInt32       stat;          // instrument status variable
	ViStatus       err;           // error number
	int            i;             // loop counter
	double         st1, st2, st3, azi, ell, dop, pow;

	// arguments
	if (argc > 1) strncpy(ip_addr, argv[1], sizeof(ip_addr));
	if (argc > 2) ip_port = (unsigned short)atoi(argv[2]);
	if (argc > 3) slot = (unsigned int)atoi(argv[3]);
	if (argc > 4) wavelength = (unsigned int)atoi(argv[4]);
	if (argc > 5) comm_tmo = (unsigned int)atoi(argv[5]);
	if (argc > 6) boot_tmo = (unsigned int)atoi(argv[6]);

	// title
	//printf("Thorlabs TXP PAX Single Measurement Sample Application\n");
	//printf("------------------------------------------------------\n\n");

	// connect to instrument
	sprintf(resName, "TCPIP::%s::%u::SOCKET", ip_addr, ip_port);   // build the resource name
	//printf("Connecting to instrument: \"%s\" Slot %u ...\n", resName, slot);
	err = TXPPAX_init(resName, slot, VI_NULL, comm_tmo, boot_tmo, 0, &instr);
	if (err)
	{
		errorMessage(VI_NULL, err);
		usage(argv[0]);
		//waitKeypress();
		return 1;
	}

	// configure measuremnt - motor speed, measurement parameters
	//printf("Configuring instrument (wavelength: %u nm) ...\n", wavelength);
	err = TXPPAX_SetOperatingMode(instr, TXPPAX_OPMODE_SINGLE_MEAS);
	if (!err) err = TXPPAX_SetWavelength(instr, (wavelength * 1.0E-9));
	if (!err) err = TXPPAX_SetPowerRange(instr, TXPPAX_RANGE_AUTO);   // autorange
	if (!err) err = TXPPAX_SetBasicMeasSpeed(instr, BASIC_MEAS_SPEED);
	if (!err) err = TXPPAX_SetSingleMeasParameters(instr, PERIODS_PER_BUFFER, BUFFERS_PER_FFT, FFT_PER_RESULT);
	if (err)
	{
		errorMessage(instr, err);
		TXPPAX_close(instr);
		//waitKeypress();
		return 1;
	}

	// check instrument status and break on status we can't measure - loop until motor has stabilized
	//printf("Checking status - waiting for motor to stabilize ...\n");
	delay1sec();
	do
	{
		err = TXPPAX_GetStatus(instr, VI_NULL, &stat);    // we don't take care of the event register value
		stat &= (TXPPAX_STATBIT_SPEED |          // ignore all other status bits
			TXPPAX_STATBIT_OVERTEMP |
			TXPPAX_STATBIT_VCC_FAIL |
			TXPPAX_STATBIT_SENSOR_MISSING |
			TXPPAX_STATBIT_SENSOR_TEMP |
			TXPPAX_STATBIT_SENSOR_ADJ |
			TXPPAX_STATBIT_CARD_ADJ);

	} while (!err && (stat == TXPPAX_STATBIT_SPEED));
	if (err)
	{
		errorMessage(instr, err);
		TXPPAX_close(instr);
		//waitKeypress();
		return 1;
	}
	if (stat)
	{
		statusMessage(stat &= ~TXPPAX_STATBIT_SPEED);	// print status message excluding motor speed failure
		TXPPAX_close(instr);
		//waitKeypress();
		return 1;
	}

	// measure in a loop and print out the results
	//printf("\n");
	//printf(" Stokes1   Stokes2   Stokes3   Azimuth    Ellipt       DOP    Power [W]\n");
	//printf("--------  --------  --------  --------  --------  --------  -----------\n");
	err = TXPPAX_RestartSingleMeas(instr);   // Clear measurement buffers
	i = 0;

	FILE *file;
	char fn[100] = "polarimeter.txt";
	file = fopen(fn, "a");
	if (file == NULL) {
		perror("Error opening file.");
	}
	char reading[1000];

	while (!err && (i < NUM_MEASUREMENTS))
	{
		err = TXPPAX_GetSingleMeas(instr, VI_TRUE, &st1, &st2, &st3, &azi, &ell, &dop, &pow, VI_NULL);
		if (!err) {
			//printf("%+8.5lf  %+8.5lf  %+8.5lf  %+8.3lf  %+8.3lf  %8.4lf  %11.3le\n", st1, st2, st3, azi, ell, dop, pow);
			sprintf(reading, "%+8.5lf,%+8.5lf,%+8.5lf,%+8.3lf,%+8.3lf,%8.4lf,%11.3le\n", st1, st2, st3, azi, ell, dop, pow);
			fprintf(file, reading);
		}
		i++;
	}
	printf(reading);

	if (err)
	{
		errorMessage(instr, err);
		TXPPAX_close(instr);
		//waitKeypress();
		return 1;
	}
	//printf("\n");

	// done - close connection
	//printf("Closing connection ...\n\n");
	TXPPAX_close(instr);
	//waitKeypress();
	return 0;
}


/*---------------------------------------------------------------------------
Error message - Query the message from the driver and print it to 'stderr'
---------------------------------------------------------------------------*/
void errorMessage(ViSession instr, ViStatus err)
{
	char buf[TXPPAX_ERR_DESCR_BUFFER_SIZE];

	TXPPAX_errorMessage(instr, err, buf);
	fprintf(stderr, "Error: %s\n\n", buf);
}


/*---------------------------------------------------------------------------
Print instrument status description - Print readable text ot status value
---------------------------------------------------------------------------*/
void statusMessage(ViUInt32 stat)
{
	printf("Status: ");
	if (stat & TXPPAX_STATBIT_CARD_EVENT)      printf("Card event occured\t");
	if (stat & TXPPAX_STATBIT_OVERTEMP)        printf("Overtemperature\t");
	if (stat & TXPPAX_STATBIT_COOLING)         printf("Cooling required\t");
	if (stat & TXPPAX_STATBIT_VCC_FAIL)        printf("Power supply failure\t");
	if (stat & TXPPAX_STATBIT_ALARM)           printf("Chassis alarm active\t");
	if (stat & TXPPAX_STATBIT_CARD_PRESENT)    printf("Card is present\t");
	if (stat & TXPPAX_STATBIT_CARD_READY)      printf("Card is ready\t");
	if (stat & TXPPAX_STATBIT_GPIO)            printf("GPIO event occured\t");
	if (stat & TXPPAX_STATBIT_SENSOR_MISSING)  printf("Sensor missing\t");
	if (stat & TXPPAX_STATBIT_MACRO_ACTIVE)    printf("Macro active\t");
	if (stat & TXPPAX_STATBIT_TRIGGER_WAIT)    printf("Waiting for trigger\t");
	if (stat & TXPPAX_STATBIT_POWER_RANGE)     printf("Input power out of range\t");
	if (stat & TXPPAX_STATBIT_SENSOR_TEMP)     printf("Sensor temperature failure\t");
	if (stat & TXPPAX_STATBIT_SPEED)           printf("Motor not stabilized\t");
	if (stat & TXPPAX_STATBIT_SENSOR_ADJ)      printf("Sensor adjustment failure\t");
	if (stat & TXPPAX_STATBIT_CARD_ADJ)        printf("Card adjustment failure\t");
	printf("\n\n");
}


/*---------------------------------------------------------------------------
Print keypress message and wait
---------------------------------------------------------------------------*/
void waitKeypress(void)
{
	printf("Press <ENTER> to exit\n");
	while (getchar() == EOF);
}


/*---------------------------------------------------------------------------
Print usage
---------------------------------------------------------------------------*/
void usage(char* name)
{
	printf("\nUsage:\n");
	printf("%s [IP-address [IP-port [slot [wavelength [comm-tmo [boot-tmo]]]]]]\n\n", name);
	printf("  IP-address    IP address or hostname (default: " DEF_IP_ADDR ")\n");
	printf("  IP-port       IP port number (default: %u)\n", DEF_IP_PORT);
	printf("  slot          zero based slot number 0..63 (default: %u)\n", DEF_SLOT);
	printf("  wavelength    wavelength [nm] (default: %u)\n", DEF_WAVELENGTH);
	printf("  comm-tmo      communication timeout [ms] (default: %u)\n", DEF_COMM_TIMEOUT);
	printf("  boot-tmo      card boot timeout [ms] (default: %u)\n", DEF_BOOT_TIMEOUT);
	printf("\n");
}


/*---------------------------------------------------------------------------
Delay ~1sec
---------------------------------------------------------------------------*/
void delay1sec(void)
{
	time_t timestamp;

	timestamp = time(NULL);
	while (time(NULL) == timestamp);
}


/****************************************************************************
End of sourc file
****************************************************************************/
