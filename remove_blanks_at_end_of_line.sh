#!/bin/sh

FF=`find . -name \*.F90 -print`
echo $FF

HF=`find . -name \*.h -print`
echo $HF

MF=`find . -name Makefile -print`
echo $MF

sed -i -e 's/[ \t]*$//' $FF $HF $MF
