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
L LM358 U?
U 1 1 599AF5BB
P 5800 4250
F 0 "U?" H 5800 4450 50  0000 L CNN
F 1 "LM358" H 5800 4050 50  0000 L CNN
F 2 "" H 5800 4250 50  0000 C CNN
F 3 "" H 5800 4250 50  0000 C CNN
	1    5800 4250
	1    0    0    -1  
$EndComp
$Comp
L D_Photo D?
U 1 1 599AF5D2
P 5000 4400
F 0 "D?" H 5020 4470 50  0000 L CNN
F 1 "D_Photo" H 4960 4290 50  0000 C CNN
F 2 "" H 4950 4400 50  0000 C CNN
F 3 "" H 4950 4400 50  0000 C CNN
	1    5000 4400
	0    1    1    0   
$EndComp
$Comp
L GND #PWR?
U 1 1 599AF718
P 5400 4650
F 0 "#PWR?" H 5400 4400 50  0001 C CNN
F 1 "GND" H 5400 4500 50  0000 C CNN
F 2 "" H 5400 4650 50  0000 C CNN
F 3 "" H 5400 4650 50  0000 C CNN
	1    5400 4650
	1    0    0    -1  
$EndComp
$Comp
L VDD #PWR?
U 1 1 599AF7C9
P 5700 3900
F 0 "#PWR?" H 5700 3750 50  0001 C CNN
F 1 "VDD" H 5700 4050 50  0000 C CNN
F 2 "" H 5700 3900 50  0000 C CNN
F 3 "" H 5700 3900 50  0000 C CNN
	1    5700 3900
	1    0    0    -1  
$EndComp
$Comp
L R R?
U 1 1 599AF851
P 5750 3400
F 0 "R?" V 5830 3400 50  0000 C CNN
F 1 "R" V 5750 3400 50  0000 C CNN
F 2 "" V 5680 3400 50  0000 C CNN
F 3 "" H 5750 3400 50  0000 C CNN
	1    5750 3400
	0    1    1    0   
$EndComp
Wire Wire Line
	6100 3400 5900 3400
Wire Wire Line
	5400 3400 5600 3400
Wire Wire Line
	6100 3150 6100 4250
Connection ~ 5400 4150
Wire Wire Line
	5400 3150 5400 4150
Wire Wire Line
	5700 3900 5700 3950
Connection ~ 5400 4650
Wire Wire Line
	5000 4500 5000 4650
Wire Wire Line
	5400 4350 5400 4650
Wire Wire Line
	5500 4350 5400 4350
Wire Wire Line
	5000 4650 5700 4650
Wire Wire Line
	5700 4650 5700 4550
Wire Wire Line
	5000 4150 5500 4150
Wire Wire Line
	5000 4200 5000 4150
$Comp
L C C?
U 1 1 599AFA6C
P 5750 3150
F 0 "C?" H 5775 3250 50  0000 L CNN
F 1 "C" H 5775 3050 50  0000 L CNN
F 2 "" H 5788 3000 50  0000 C CNN
F 3 "" H 5750 3150 50  0000 C CNN
	1    5750 3150
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5900 3150 6100 3150
Connection ~ 6100 3400
Wire Wire Line
	5600 3150 5400 3150
Connection ~ 5400 3400
Wire Wire Line
	6100 4250 6450 4250
$EndSCHEMATC
