#!/bin/bash

VAR=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/variant.py
SIFT=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/sift.py
SEATTLE=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/seattle.py
SNP=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/snp.py
ANN=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/annotator.py
EFF=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/snpeff.py
BROAD=aheiberg@triton-login.sdsc.edu:/home/aheiberg/andrew/broad.py
rsync -avP $VAR $SIFT $SEATTLE $SNP $ANN $EFF $BROAD .
