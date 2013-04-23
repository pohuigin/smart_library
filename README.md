smart_library
=============

An SSW IDL routine library to be used by SMART itself and users wishing to analyse SMART detections.

Requirements:

Add the following to your IDL start-up file (or run manually before using SMART_LIBRARY routines):

;Set up system variables for SMART_Library codes

;------------------------------------------------------------------------------>

;Set the environment variables

DEFSYSV, '!AR_PATH', '~/science/repositories/smart_library/'

DEFSYSV, '!AR_PARAM', 'ar_param.txt'

;------------------------------------------------------------------------------>

Notes:

1. This repository requires an SSWIDL installation

2. This repository has dependencies as listed below

Dependencies:

	GEN_LIBRARY: git@github.com:pohuigin/gen_library.git
