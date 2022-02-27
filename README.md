# MUS-177-HW2 by Angus Yick

This is the second assignment for my music programming class. 

![image](https://user-images.githubusercontent.com/74380180/155871711-d4203bae-241a-47d1-829d-15cf3a348a72.png)

```
wavefold_allpass~.dll: wavefold_allpass~.c; 
	cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	link /dll /export:wavefold_allpass_tilde_setup $*.obj $(PDNTLIB)
```
