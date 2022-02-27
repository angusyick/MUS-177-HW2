# MUS-177-HW2 by Angus Yick

This is the second assignment for my music programming class. It is a PD external that shows some techniques such as: lookup table oscillator, wavefolding, and allpass filter. Functionality of the external is outlined in the help patch. 

![image](https://user-images.githubusercontent.com/74380180/155871711-d4203bae-241a-47d1-829d-15cf3a348a72.png)

Add the following to the makefile:

In pd_nt:
```
wavefold_allpass~.dll
```

At the bottom:
```
wavefold_allpass~.dll: wavefold_allpass~.c; 
	cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	link /dll /export:wavefold_allpass_tilde_setup $*.obj $(PDNTLIB)
```
