
// exec builder.sce  //  must be run from this directory to create mex function library
// It is expected that ../../mex and ../../src  are copies of the corresponding
// OpenFEM_Matlab directories
//
//
//
// To save directory 
// zip -r of_scilabcvs * -x */CVS/* */*/CVS/* openfem_scilab/mex/*.*o */*.*o src/*.f openfem_scilab/.libs/*
// cd ~/of_scilab; rm of_scilabcvs.zip
//
// Etienne Balmes, Anne-Sophie Mouronval
// $Revision: 1.1 $  $Date: 2006/01/18 13:19:16 $ 


// cleaning links and other
lines(0);
unix("ln -s ../src/*.c .");
unix("rm [a-k]*.o [a-k]*.lo [m-z]*.o [m-z]*.lo lib_of*");


ilib_name  = 'lib_ofmex' 		// interface library name 
files = ['tab0d.o';'tab1d.o';'tab2d.o';'taba0d.o';'taba1d.o';'taba2d.o';'taba6d.o';'taba8d.o';'tabaxd.o';'tabszd.o';'ab0d.o';'ab1d.o';'ab4d.o';'ab5d.o';'hookax.o';'plmasd.o';'plonad.o';'etr3p1d.o';'etm3p1d.o';'ets3p1d.o';'etc3p1d.o';'etr3q1d.o';'etm3q1d.o';'ets3q1d.o';'etc3q1d.o';'etm3q2c.o';'etr3q2c.o';'ets3q2c.o';'ec3c2c.o';'etc3q2c.o';'etr3r1d.o';'etm3r1d.o';'ets3r1d.o';'etc3r1d.o';'er3c2c.o';'em3c2c.o';'es3d2c.o';'es3c2c.o';'emaq2c.o';'dpaq2c.o';'etr3r2c.o';'etm3r2c.o';'ets3r2c.o';'etc3r2c.o';'etr3p2c.o';'etm3p2c.o';'ets3p2c.o';'etc3p2c.o';'fobase.o';'dcopy.o';'etm2p1d.o';'etr2p1d.o';'ets2p1d.o';'etc2p1d.o';'etm2p2c.o';'etr2p2c.o';'ets2p2c.o';'etc2p2c.o';'etm2q1d.o';'etr2q1d.o';'ets2q1d.o';'etc2q1d.o';'etm2q2c.o';'etr2q2c.o';'ets2q2c.o';'etc2q2c.o';'e1ap1d.o';'etrap1d.o';'etmap1d.o';'etsap1d.o';'etcap1d.o';'e2ap2c.o';'etrap2c.o';'etmap2c.o';'etsap2c.o';'etcap2c.o';'e1aq1c.o';'etmaq1d.o';'etraq1d.o';'etsaq1d.o';'etcaq1d.o';'e2aq2c.o';'etmaq2c.o';'etraq2c.o';'eraq2c.o';'etsaq2c.o';'etcaq2c.o';'etmdktp.o';'etrdktp.o';'etsdktp.o';'etcdktp.o';'etm5noe.o';'etr5noe.o';'ets5noe.o';'etc5noe.o';'etrmit4.o';'etsmit4.o';'bremit.o';'etsmitx.o';'inmit4.o';'mitlin.o';'fonmit.o';'melmit.o';'replo1.o';'pcq12d.o';'canoq1.o';'chan56.o';'ddot.o';'of_mk_subs.o';'chan57mod.o'];    // objects files 
				// 
libs  = ['libscs'] 				// other libs needed for linking
makename = ['Makelib']
ldflags = []

// xxx missing an automated selection of OSTYPE xxx 

cflags = ['-DOSTYPEmexglx']                // C compilation options
fflags = []
table = [ 'of_mk', 'of_mk' 'cmex';
          'of_time' 'of_time' 'cmex';
          'sp_util' 'sp_util' 'cmex';
          'nopo2sd' 'nopo2sd' 'cmex']  // table of (scilab_name,interface-name) 

// ----------------------------------------------
// now build the library
ilib_build(ilib_name,table,files,libs,makename,ldflags,cflags,fflags)

// ----------------------------------------------
// Finally load it for test purposes
exec loader.sce
of_mk cvs
of_time cvs
nopo2sd cvs
sp_util cvs

