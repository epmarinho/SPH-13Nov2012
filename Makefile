#!/bin/make

CFLAGS=-I./include

target=default
$(target):
	(cd src;make $@)