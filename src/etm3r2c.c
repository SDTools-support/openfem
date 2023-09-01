/* etm3r2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__15 = 15;
static int32 c__21 = 21;

/* Subroutine */ int etm3r2c_(coor, ro, iopt, ae)
doublereal *coor, *ro;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 ijt[45] = { 1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45 };
    static doublereal poids[21] = { .0174915528534482,.0174915528534482,.0174915528534482,.0183880767761814,.0183880767761814,.0183880767761814,.0312499991721577,.0279864845655171,.0279864845655171,.0279864845655171,.0294209228418903,.0294209228418903,.0294209228418903,.0499999986754523,.0174915528534482,.0174915528534482,.0174915528534482,.0183880767761814,.0183880767761814,.0183880767761814,.0312499991721577 };
    static struct {
	doublereal e_1[945];
	} equiv_6 = { -1.74292420972099, -1.74292420972099, -1.7097211822211, -.727812937550124, 0., -.0761437882318176, 0., -.727812937550124, -.0761437882318176, -.0467837316913592, -.0467837316913592, -.761015965050019, -.26704103315605, 0., -.237680976615592, 0., -.26704103315605, -.237680976615592, 2.47073714727111, -.359485397070617, -.323074376767549, .359485397070617, .359485397070617, -.0410358262631383, -.359485397070617, 2.47073714727111, -.323074376767549, -.4, -.4, 2.47073714727111, 
		.4, 0., .31382476484741, 0., .4, .31382476484741, .31382476484741, -.0456606322232079, .323074376767549, .0456606322232079, .0456606322232079, .0410358262631383, -.0456606322232079, .31382476484741, .323074376767549, .727812937550124, .727812937550124, -.0761437882318179, 1.74292420972099, 0., -1.7097211822211, 0., -.727812937550124, -.0761437882318176, .26704103315605, .26704103315605, -.237680976615592, .0467837316913591, 0., -.761015965050018, 0., -.26704103315605, 
		-.237680976615592, -2.47073714727111, -2.83022254434173, -.323074376767549, .359485397070617, 2.83022254434173, -.323074376767549, -.359485397070617, 0., -.0410358262631383, -.4, -.4, .31382476484741, .4, 0., 2.47073714727111, 0., .4, .31382476484741, -.313824764847409, -.359485397070617, .323074376767549, .0456606322232079, .359485397070617, .323074376767549, -.0456606322232079, 0., .0410358262631383, .727812937550124, .727812937550124, -.0761437882318177, -.727812937550124, 0., 
		-.0761437882318176, 0., 1.74292420972099, -1.7097211822211, .26704103315605, .26704103315605, -.237680976615592, -.26704103315605, 0., -.237680976615592, 0., .0467837316913591, -.761015965050018, 0., -.359485397070617, -.0410358262631383, 2.83022254434173, .359485397070617, -.323074376767549, -2.83022254434173, -2.47073714727111, -.323074376767549, -.4, -.4, .31382476484741, .4, 0., .31382476484741, 0., .4, 2.47073714727111, 0., -.0456606322232079, .0410358262631383, 
		.359485397070617, .0456606322232079, .323074376767549, -.359485397070617, -.313824764847409, .323074376767549, -.581326747441765, -.581326747441765, -.700266010629198, .581326747441765, 0., -.700266010629197, 0., -.875355160262788, -.0399275296758689, .100758491021304, .100758491021304, -.756415897075355, -.100758491021304, 0., -.756415897075355, 0., -.285781352578133, -.14509533188096, 0., -1.66862508206251, -.884134241764073, .211943174357954, 1.66862508206251, 
		-.112299772892315, -.211943174357954, 1.45668190770455, -.112299772892315, -.4, -.4, 1.45668190770455, .4, 0., 1.45668190770455, 0., .4, .185022861556829, 0., -.211943174357954, .884134241764073, .0269203128011253, .211943174357954, .112299772892315, -.0269203128011253, .185022861556829, .112299772892315, .875355160262787, .875355160262787, -.0399275296758689, .581326747441765, 0., -.700266010629197, 0., .581326747441765, -.700266010629197, .285781352578133, .285781352578133, 
		-.14509533188096, -.100758491021304, 0., -.756415897075355, 0., -.100758491021304, -.756415897075355, -1.45668190770455, -1.66862508206251, -.112299772892315, 1.66862508206251, 1.66862508206251, -.884134241764072, -1.66862508206251, -1.45668190770455, -.112299772892315, -.4, -.4, .185022861556829, .4, 0., 1.45668190770455, 0., .4, 1.45668190770455, -.185022861556829, -.211943174357954, .112299772892315, .211943174357954, .211943174357954, .884134241764072, -.211943174357954, 
		-.185022861556829, .112299772892315, -.581326747441765, -.581326747441765, -.700266010629198, -.875355160262788, 0., -.0399275296758689, 0., .581326747441765, -.700266010629197, .100758491021304, .100758491021304, -.756415897075355, -.285781352578133, 0., -.14509533188096, 0., -.100758491021304, -.756415897075355, 1.45668190770455, -.211943174357954, -.112299772892315, 1.66862508206251, .211943174357954, -.112299772892315, -1.66862508206251, 0., -.884134241764073, -.4, -.4, 
		1.45668190770455, .4, 0., .185022861556829, 0., .4, 1.45668190770455, .185022861556829, -.0269203128011253, .112299772892315, .211943174357954, .0269203128011253, .112299772892315, -.211943174357954, 0., .884134241764073, -.095766041024111, -.095766041024111, -.405286630980767, .0957661467983153, 0., -.405286687084434, 0., .0957661467983153, -.405286687084434, .162432787163638, .162432787163638, -.627508866448465, -.162432773728552, 0., -.627508902683917, 0., -.162432773728552, 
		-.627508902683917, -1.05774204062214e-7, -1.18306448141906, -.444444431198967, 1.18306448141906, 1.18306448141906, -.444444470935398, -1.18306448141906, -1.05774204062214e-7, -.444444431198967, -.4, -.4, 1.03279549742923, .4, 0., 1.03279558976835, 0., .4, 1.03279558976835, -1.34350854885668e-8, -.150268891650706, .444444431198967, .150268891650706, .150268891650706, .444444470935398, -.150268891650706, -1.34350854885668e-8, .444444431198967, -.594853970706175, -.594853970706175, 
		-.474352608585539, -.797426985353087, 0., .0807685941918872, 0., -.797426985353087, .0807685941918872, -.594853970706175, -.594853970706175, .474352608585539, -.797426985353087, 0., -.0807685941918872, 0., -.797426985353087, -.0807685941918872, 1.39228095605926, -.202573014646913, -.323074376767549, .202573014646913, .202573014646913, -.0410358262631383, -.202573014646913, 1.39228095605926, -.323074376767549, -1., -1., 0., 1., 0., 0., 0., 1., 0., 1.39228095605926, 
		-.202573014646913, .323074376767549, .202573014646913, .202573014646913, .0410358262631383, -.202573014646913, 1.39228095605926, .323074376767549, .797426985353087, .797426985353087, .0807685941918874, .594853970706174, 0., -.474352608585538, 0., -.797426985353087, .0807685941918872, .797426985353087, .797426985353087, -.080768594191887, .594853970706174, 0., .474352608585539, 0., -.797426985353087, -.0807685941918872, -1.39228095605926, -1.59485397070617, -.323074376767549, 
		.202573014646913, 1.59485397070617, -.323074376767549, -.202573014646913, 0., -.0410358262631383, -1., -1., 0., 1., 0., 0., 0., 1., 0., -1.39228095605926, -1.59485397070617, .323074376767549, .202573014646913, 1.59485397070617, .323074376767549, -.202573014646913, 0., .0410358262631383, .797426985353087, .797426985353087, .0807685941918872, -.797426985353087, 0., .0807685941918872, 0., .594853970706174, -.474352608585538, .797426985353087, .797426985353087, -.0807685941918871, 
		-.797426985353087, 0., -.0807685941918872, 0., .594853970706174, .474352608585539, 0., -.202573014646913, -.0410358262631383, 1.59485397070617, .202573014646913, -.323074376767549, -1.59485397070617, -1.39228095605926, -.323074376767549, -1., -1., 0., 1., 0., 0., 0., 1., 0., 0., -.202573014646913, .0410358262631383, 1.59485397070617, .202573014646913, .323074376767549, -1.59485397070617, -1.39228095605926, .323074376767549, .0597158717897698, .0597158717897698, .0280749432230789, 
		-.0597158717897699, 0., .0280749432230789, 0., -.88056825642046, .0525839011025453, .0597158717897697, .0597158717897697, -.0280749432230787, -.0597158717897699, 0., -.0280749432230789, 0., -.88056825642046, -.0525839011025453, 0., -.94028412821023, -.884134241764073, .11943174357954, .94028412821023, -.112299772892315, -.11943174357954, .820852384630691, -.112299772892315, -1., -1., 0., 1., 0., 0., 0., 1., 0., 0., -.94028412821023, .884134241764073, .11943174357954, 
		.94028412821023, .112299772892315, -.11943174357954, .820852384630691, .112299772892315, .88056825642046, .88056825642046, .0525839011025456, -.0597158717897699, 0., .0280749432230789, 0., -.0597158717897699, .0280749432230789, .88056825642046, .88056825642046, -.0525839011025456, -.0597158717897699, 0., -.0280749432230789, 0., -.0597158717897699, -.0280749432230789, -.82085238463069, -.94028412821023, -.112299772892315, .94028412821023, .94028412821023, -.884134241764072, 
		-.94028412821023, -.82085238463069, -.112299772892315, -1., -1., 0., 1., 0., 0., 0., 1., 0., -.82085238463069, -.94028412821023, .112299772892315, .94028412821023, .94028412821023, .884134241764072, -.94028412821023, -.82085238463069, .112299772892315, .0597158717897697, .0597158717897697, .0280749432230788, -.88056825642046, 0., .0525839011025453, 0., -.0597158717897699, .0280749432230789, .0597158717897697, .0597158717897697, -.0280749432230788, -.88056825642046, 0., 
		-.0525839011025453, 0., -.0597158717897699, -.0280749432230789, .820852384630691, -.11943174357954, -.112299772892315, .94028412821023, .11943174357954, -.112299772892315, -.94028412821023, 0., -.884134241764073, -1., -1., 0., 1., 0., 0., 0., 1., 0., .820852384630691, -.11943174357954, .112299772892315, .94028412821023, .11943174357954, .112299772892315, -.94028412821023, 0., .884134241764073, .333333373069763, .333333373069763, .111111117733849, -.333333313465118, 0., 
		.111111107799742, 0., -.333333313465118, .111111107799742, .333333373069763, .333333373069763, -.111111117733849, -.333333313465118, 0., -.111111107799742, 0., -.333333313465118, -.111111107799742, -5.96046447753906e-8, -.666666686534882, -.444444431198967, .666666686534882, .666666686534882, -.444444470935398, -.666666686534882, -5.96046447753906e-8, -.444444431198967, -1., -1., 0., 1., 0., 0., 0., 1., 0., -5.96046447753906e-8, -.666666686534882, .444444431198967, 
		.666666686534882, .666666686534882, .444444470935398, -.666666686534882, -5.96046447753906e-8, .444444431198967, -.0467837316913596, -.0467837316913596, .761015965050018, -.26704103315605, 0., .237680976615592, 0., -.26704103315605, .237680976615592, -1.74292420972099, -1.74292420972099, 1.7097211822211, -.727812937550124, 0., .0761437882318176, 0., -.727812937550124, .0761437882318176, .31382476484741, -.0456606322232079, -.323074376767549, .0456606322232079, .0456606322232079, 
		-.0410358262631383, -.0456606322232079, .31382476484741, -.323074376767549, -.4, -.4, -2.47073714727111, .4, 0., -.31382476484741, 0., .4, -.31382476484741, 2.47073714727111, -.359485397070617, .323074376767549, .359485397070617, .359485397070617, .0410358262631383, -.359485397070617, 2.47073714727111, .323074376767549, .26704103315605, .26704103315605, .237680976615592, .0467837316913591, 0., .761015965050018, 0., -.26704103315605, .237680976615592, .727812937550123, 
		.727812937550123, .0761437882318183, 1.74292420972099, 0., 1.7097211822211, 0., -.727812937550124, .0761437882318176, -.31382476484741, -.359485397070618, -.323074376767549, .0456606322232079, .359485397070618, -.323074376767549, -.0456606322232079, 0., -.0410358262631383, -.4, -.4, -.31382476484741, .4, 0., -2.47073714727111, 0., .4, -.31382476484741, -2.47073714727111, -2.83022254434173, .323074376767549, .359485397070617, 2.83022254434173, .323074376767549, -.359485397070617, 
		0., .0410358262631383, .26704103315605, .26704103315605, .237680976615592, -.26704103315605, 0., .237680976615592, 0., .0467837316913591, .761015965050018, .727812937550123, .727812937550123, .0761437882318179, -.727812937550124, 0., .0761437882318176, 0., 1.74292420972099, 1.7097211822211, 0., -.0456606322232079, -.0410358262631383, .359485397070618, .0456606322232079, -.323074376767549, -.359485397070618, -.31382476484741, -.323074376767549, -.4, -.4, -.313824764847411, .4, 0., 
		-.31382476484741, 0., .4, -2.47073714727111, 0., -.359485397070617, .0410358262631383, 2.83022254434173, .359485397070617, .323074376767549, -2.83022254434173, -2.47073714727111, .323074376767549, .100758491021304, .100758491021304, .756415897075355, -.100758491021304, 0., .756415897075355, 0., -.285781352578133, .14509533188096, -.581326747441765, -.581326747441765, .700266010629197, .581326747441765, 0., .700266010629197, 0., -.875355160262788, .039927529675869, 0., 
		-.211943174357954, -.884134241764073, .0269203128011253, .211943174357954, -.112299772892315, -.0269203128011253, .185022861556829, -.112299772892315, -.4, -.4, -1.45668190770455, .4, 0., -1.45668190770455, 0., .4, -.185022861556829, 0., -1.66862508206251, .884134241764073, .211943174357954, 1.66862508206251, .112299772892315, -.211943174357954, 1.45668190770455, .112299772892315, .285781352578132, .285781352578132, .14509533188096, -.100758491021304, 0., .756415897075355, 0., 
		-.100758491021304, .756415897075355, .875355160262787, .875355160262787, .0399275296758687, .581326747441765, 0., .700266010629197, 0., .581326747441765, .700266010629197, -.185022861556828, -.211943174357954, -.112299772892315, .211943174357954, .211943174357954, -.884134241764072, -.211943174357954, -.185022861556828, -.112299772892315, -.4, -.4, -.185022861556829, .4, 0., -1.45668190770455, 0., .4, -1.45668190770455, -1.45668190770455, -1.66862508206251, .112299772892315, 
		1.66862508206251, 1.66862508206251, .884134241764072, -1.66862508206251, -1.45668190770455, .112299772892315, .100758491021304, .100758491021304, .756415897075355, -.285781352578133, 0., .14509533188096, 0., -.100758491021304, .756415897075355, -.581326747441765, -.581326747441765, .700266010629197, -.875355160262788, 0., .039927529675869, 0., .581326747441765, .700266010629197, .185022861556829, -.0269203128011253, -.112299772892315, .211943174357954, .0269203128011253, 
		-.112299772892315, -.211943174357954, 0., -.884134241764073, -.4, -.4, -1.45668190770455, .4, 0., -.185022861556829, 0., .4, -1.45668190770455, 1.45668190770455, -.211943174357954, .112299772892315, 1.66862508206251, .211943174357954, .112299772892315, -1.66862508206251, 0., .884134241764073, .162432787163637, .162432787163637, .627508866448465, -.162432773728552, 0., .627508902683917, 0., -.162432773728552, .627508902683917, -.0957660410241115, -.0957660410241115, 
		.405286630980767, .0957661467983153, 0., .405286687084434, 0., .0957661467983153, .405286687084434, -1.34350854885668e-8, -.150268891650706, -.444444431198967, .150268891650706, .150268891650706, -.444444470935398, -.150268891650706, -1.34350854885668e-8, -.444444431198967, -.4, -.4, -1.03279549742923, .4, 0., -1.03279558976835, 0., .4, -1.03279558976835, -1.05774204062214e-7, -1.18306448141906, .444444431198967, 1.18306448141906, 1.18306448141906, .444444470935398, 
		-1.18306448141906, -1.05774204062214e-7, .444444431198967 };

    static struct {
	doublereal e_1[315];
	} equiv_9 = { .261406882550335, -.0919231405808113, -.0919231405808113, -.106025068106032, -.0293600565404584, -.0293600565404584, .28666335646448, .0364110203030687, .28666335646448, .318970794141235, .0405146029293825, .0405146029293825, .0364110203030687, .00462480596006959, .0364110203030687, -.0919231405808111, .261406882550335, -.0919231405808113, -.0293600565404585, -.106025068106032, -.0293600565404584, .28666335646448, .28666335646448, .0364110203030687, .0405146029293826, 
		.318970794141235, .0405146029293825, .0364110203030688, .0364110203030687, .0046248059600696, -.0919231405808111, -.0919231405808113, .261406882550335, -.0293600565404585, -.0293600565404584, -.106025068106032, .0364110203030687, .28666335646448, .28666335646448, .0405146029293826, .0405146029293825, .318970794141235, .0046248059600696, .0364110203030687, .0364110203030688, -.118939263187433, -.118939263187433, -.0586007822341042, -.0971925056776921, -.0971925056776921, 
		-.017869467584349, .784490840298434, .0996434014656387, .0996434014656388, .188056825642046, .188056825642046, .0238863487159079, .0996434014656388, .0126563714266764, .0126563714266764, -.0586007822341042, -.118939263187433, -.118939263187433, -.0178694675843491, -.0971925056776921, -.0971925056776921, .0996434014656391, .784490840298434, .0996434014656391, .0238863487159079, .188056825642046, .188056825642046, .0126563714266764, .0996434014656387, .0126563714266764, 
		-.118939263187433, -.0586007822341042, -.118939263187433, -.0971925056776921, -.017869467584349, -.0971925056776921, .0996434014656388, .0996434014656387, .784490840298434, .188056825642046, .0238863487159079, .188056825642046, .0126563714266764, .0126563714266764, .0996434014656388, -.165255372416117, -.165255369562065, -.165255369562065, -.0791890707037792, -.0791890755446534, -.0791890755446534, .394354803634306, .394354838892375, .394354803634306, .133333325386047, 
		.133333337306976, .133333337306976, .0500896275646608, .0500896320430227, .0500896275646608, -.161537188383774, -.0910275507576717, -.0910275507576717, -.161537188383774, -.0910275507576717, -.0910275507576717, .161537188383774, .0205179131315691, .161537188383774, .797426985353087, .101286507323456, .101286507323456, .161537188383774, .0205179131315691, .161537188383774, -.0910275507576718, -.161537188383774, -.0910275507576717, -.0910275507576718, -.161537188383774, 
		-.0910275507576717, .161537188383775, .161537188383774, .0205179131315692, .101286507323456, .797426985353087, .101286507323456, .161537188383774, .161537188383774, .0205179131315692, -.0910275507576717, -.0910275507576717, -.161537188383774, -.0910275507576718, -.0910275507576717, -.161537188383774, .0205179131315692, .161537188383774, .161537188383774, .101286507323456, .101286507323456, .797426985353087, .0205179131315692, .161537188383774, .161537188383774, -.249108503664097, 
		-.249108503664097, -.0561498864461576, -.249108503664097, -.249108503664097, -.0561498864461576, .442067120882036, .0561498864461576, .0561498864461576, .470142064105115, .470142064105115, .0597158717897698, .442067120882036, .0561498864461576, .0561498864461576, -.0561498864461577, -.249108503664097, -.249108503664097, -.0561498864461577, -.249108503664097, -.249108503664097, .0561498864461577, .442067120882036, .0561498864461577, .0597158717897699, .470142064105115, 
		.470142064105115, .0561498864461577, .442067120882036, .0561498864461577, -.249108503664097, -.0561498864461576, -.249108503664097, -.249108503664097, -.0561498864461576, -.249108503664097, .0561498864461576, .0561498864461576, .442067120882036, .470142064105115, .0597158717897698, .470142064105115, .0561498864461576, .0561498864461576, .442067120882036, -.222222215599484, -.222222225533591, -.222222225533591, -.222222215599484, -.222222225533591, -.222222225533591, 
		.222222215599484, .222222235467699, .222222215599484, .333333313465118, .333333343267441, .333333343267441, .222222215599484, .222222235467699, .222222215599484, -.106025068106032, -.0293600565404584, -.0293600565404584, .261406882550335, -.0919231405808112, -.0919231405808112, .0364110203030687, .00462480596006959, .0364110203030687, .318970794141235, .0405146029293825, .0405146029293825, .28666335646448, .0364110203030687, .28666335646448, -.0293600565404584, -.106025068106032, 
		-.0293600565404584, -.0919231405808112, .261406882550336, -.0919231405808112, .0364110203030688, .0364110203030687, .00462480596006959, .0405146029293828, .318970794141235, .0405146029293825, .28666335646448, .28666335646448, .0364110203030688, -.0293600565404582, -.0293600565404584, -.106025068106032, -.0919231405808113, -.0919231405808112, .261406882550336, .00462480596006959, .0364110203030687, .0364110203030688, .0405146029293826, .0405146029293825, .318970794141235, 
		.0364110203030687, .28666335646448, .28666335646448, -.0971925056776921, -.0971925056776922, -.017869467584349, -.118939263187433, -.118939263187433, -.0586007822341042, .0996434014656387, .0126563714266764, .0126563714266764, .188056825642046, .188056825642046, .0238863487159079, .784490840298434, .0996434014656387, .0996434014656387, -.0178694675843491, -.0971925056776922, -.0971925056776922, -.0586007822341045, -.118939263187433, -.118939263187433, .0126563714266763, 
		.0996434014656388, .0126563714266763, .0238863487159078, .188056825642046, .188056825642046, .0996434014656391, .784490840298434, .0996434014656391, -.0971925056776922, -.017869467584349, -.0971925056776922, -.118939263187433, -.0586007822341042, -.118939263187433, .0126563714266764, .0126563714266764, .0996434014656387, .188056825642046, .0238863487159079, .188056825642046, .0996434014656388, .0996434014656387, .784490840298434, -.079189070703779, -.0791890755446533, 
		-.0791890755446533, -.165255372416117, -.165255369562065, -.165255369562065, .0500896275646607, .0500896320430227, .0500896275646607, .133333325386047, .133333337306976, .133333337306976, .394354803634306, .394354838892375, .394354803634306 };


    /* Local variables */
    extern /* Subroutine */ int em3c2c_();
    static doublereal delta[21];
#define dp1 ((doublereal *)&equiv_6)
#define dp2 ((doublereal *)&equiv_6 + 135)
#define dp3 ((doublereal *)&equiv_6 + 270)
#define dp4 ((doublereal *)&equiv_6 + 405)
#define dp5 ((doublereal *)&equiv_6 + 540)
#define dp6 ((doublereal *)&equiv_6 + 675)
#define dp7 ((doublereal *)&equiv_6 + 810)
#define va1 ((doublereal *)&equiv_9)
#define va2 ((doublereal *)&equiv_9 + 105)
#define va3 ((doublereal *)&equiv_9 + 210)
#define vp1 ((doublereal *)&equiv_9)
#define dpr ((doublereal *)&equiv_6)

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE masse */
/*  ---     pentaedre r1 DROIT */
/*          formule d'integration a 21 points (gauss) 7*3 */

/*  PARAMETRES D ENTREE  : */
/*  ------------------- */
/*   X,Y,Z   : TABLEAUX DES COORDONNEES DES POINTS DE L ELEMENT */
/*   ijt     : permutation pour oasser de la numerotation par inconnues */
/*             a celle par noeuds */
/*   nno     : nombre de noeuds de l'element */
/*   npo     : nombre de points */
/*   npi     : nombre de points d'integration */
/*   dp      : valeur des derivees des polynomes de base aux points */
/*             d'integration sur l'element de reference */
/*   vp1     : valeur des polynomes de base aux points */
/*             d'integration sur l'element de reference */
/*   poids   : poids de la formule d'integration */

/*  tableaux de travail : */
/*  ------------------- */
/*   delta   : jacobien aux points d'integration */
/*   (x y z)int : coordonnees des points d'integration sur */
/*              l'element courant */

/*  PARAMETRE DE SORTIE  : */
/*  -------------------- */
/*   AE     : MATRICE DE Masse . SYMETRIQUE */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    --ae;
    coor -= 16;

    /* Function Body */


    em3c2c_(&c__15, &c__15, &coor[16], &coor[31], &coor[46], &c__21, ijt, poids, vp1, dpr, ro, &ae[1], delta);
} /* etm3r2c_ */

#undef dpr
#undef vp1
#undef va3
#undef va2
#undef va1
#undef dp7
#undef dp6
#undef dp5
#undef dp4
#undef dp3
#undef dp2
#undef dp1


