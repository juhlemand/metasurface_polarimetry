#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "TXP_DRV_PAX.h"


void errorMessage(ViSession instr, ViStatus err)
{
	char buf[TXPPAX_ERR_DESCR_BUFFER_SIZE];

	TXPPAX_errorMessage(instr, err, buf);
	fprintf(stderr, "Error: %s\n\n", buf);
}

int main(int argc, char** argv) {
	char           resName[80];
	char           ip_addr[60] = "localhost"; //localhost
	unsigned short ip_port = 2402;
	unsigned short slotNumber = 0;
	unsigned int wavelength = 532;
	double         st1, st2, st3, azi, ell, dop, pow;

	ViSession      instr;         // instruenmt handle
	ViUInt32       stat;          // instrument status variable
	ViStatus       err;           // error number
	ViUInt32 communicationTimeout = 60000;
	ViUInt32 cardBootTimeout = 60000;
	ViUInt32 keepAliveTime = 0;

	sprintf(resName, "TCPIP::%s::%u::SOCKET", ip_addr, ip_port);
	printf(resName);
	err = TXPPAX_init(resName, slotNumber, VI_NULL, communicationTimeout, cardBootTimeout, keepAliveTime, &instr);
	if (err)
	{
		errorMessage(VI_NULL, err);
		TXPPAX_close(instr);
		return 1;
	}
	TXPPAX_SetOperatingMode(instr, TXPPAX_OPMODE_SINGLE_MEAS);
	TXPPAX_SetWavelength(instr, (wavelength * 1.0E-9));
	//TXPPAX_SetBasicMeasSpeed(instr, 200.0);
	TXPPAX_GetStatus(instr, VI_NULL, &stat);


	FILE *file;
	char fn[100] = "polarimeter.txt";

	//if (file == NULL) {
	//	perror("Error opening file.");
	//}
	int i;
	char reading[1000];
	Sleep(500);
	while(1==1) {
		TXPPAX_GetSingleMeas(instr, VI_TRUE, &st1, &st2, &st3, &azi, &ell, &dop, &pow, VI_NULL);
		file = fopen(fn, "a");
		sprintf(reading, "%+8.5lf,%+8.5lf,%+8.5lf,%+8.3lf,%+8.3lf,%8.4lf,%11.3le\n", st1, st2, st3, azi, ell, dop, pow);
		//printf(reading);
		fprintf(file, reading);
		fclose(file);
		Sleep(100);
	}
	//printf(reading);
	

	Sleep(500);
	TXPPAX_close(instr);
	return 0;
}
