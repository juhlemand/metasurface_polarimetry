EESchema Schematic File Version 2
LIBS:power
LIBS:device
LIBS:transistors
LIBS:conn
LIBS:linear
LIBS:regul
LIBS:74xx
LIBS:cmos4000
LIBS:adc-dac
LIBS:memory
LIBS:xilinx
LIBS:microcontrollers
LIBS:dsp
LIBS:microchip
LIBS:analog_switches
LIBS:motorola
LIBS:texas
LIBS:intel
LIBS:audio
LIBS:interface
LIBS:digital-audio
LIBS:philips
LIBS:display
LIBS:cypress
LIBS:siliconi
LIBS:opto
LIBS:atmel
LIBS:contrib
LIBS:valves
EELAYER 25 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L LM358 U1
U 1 1 599AF5BB
P 5400 4100
F 0 "U1" H 5400 4300 50  0000 L CNN
F 1 "LM358" H 5400 3900 50  0000 L CNN
F 2 "" H 5400 4100 50  0000 C CNN
F 3 "" H 5400 4100 50  0000 C CNN
	1    5400 4100
	1    0    0    -1  
$EndComp
$Comp
L D_Photo D1
U 1 1 599AF5D2
P 4600 4250
F 0 "D1" H 4620 4320 50  0000 L CNN
F 1 "D_Photo" H 4560 4140 50  0000 C CNN
F 2 "" H 4550 4250 50  0000 C CNN
F 3 "" H 4550 4250 50  0000 C CNN
	1    4600 4250
	0    1    1    0   
$EndComp
$Comp
L GND #PWR1
U 1 1 599AF718
P 5000 4500
F 0 "#PWR1" H 5000 4250 50  0001 C CNN
F 1 "GND" H 5000 4350 50  0000 C CNN
F 2 "" H 5000 4500 50  0000 C CNN
F 3 "" H 5000 4500 50  0000 C CNN
	1    5000 4500
	1    0    0    -1  
$EndComp
$Comp
L VDD #PWR2
U 1 1 599AF7C9
P 5300 3750
F 0 "#PWR2" H 5300 3600 50  0001 C CNN
F 1 "VDD" H 5300 3900 50  0000 C CNN
F 2 "" H 5300 3750 50  0000 C CNN
F 3 "" H 5300 3750 50  0000 C CNN
	1    5300 3750
	1    0    0    -1  
$EndComp
$Comp
L R R1
U 1 1 599AF851
P 5350 3250
F 0 "R1" V 5430 3250 50  0000 C CNN
F 1 "R" V 5350 3250 50  0000 C CNN
F 2 "" V 5280 3250 50  0000 C CNN
F 3 "" H 5350 3250 50  0000 C CNN
	1    5350 3250
	0    1    1    0   
$EndComp
Wire Wire Line
	5700 3250 5500 3250
Wire Wire Line
	5000 3250 5200 3250
Wire Wire Line
	5700 3000 5700 4100
Connection ~ 5000 4000
Wire Wire Line
	5000 3000 5000 4000
Wire Wire Line
	5300 3750 5300 3800
Connection ~ 5000 4500
Wire Wire Line
	4600 4350 4600 4500
Wire Wire Line
	5000 4200 5000 4500
Wire Wire Line
	5100 4200 5000 4200
Wire Wire Line
	4600 4500 5300 4500
Wire Wire Line
	5300 4500 5300 4400
Wire Wire Line
	4600 4000 5100 4000
Wire Wire Line
	4600 4050 4600 4000
$Comp
L C C1
U 1 1 599AFA6C
P 5350 3000
F 0 "C1" H 5375 3100 50  0000 L CNN
F 1 "C" H 5375 2900 50  0000 L CNN
F 2 "" H 5388 2850 50  0000 C CNN
F 3 "" H 5350 3000 50  0000 C CNN
	1    5350 3000
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5500 3000 5700 3000
Connection ~ 5700 3250
Wire Wire Line
	5200 3000 5000 3000
Connection ~ 5000 3250
Wire Wire Line
	5700 4100 6050 4100
Text Label 6050 4100 0    60   ~ 0
ADC
$EndSCHEMATC
