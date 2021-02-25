# blacktern

Black tern is a software used to calculate the wave velocity values of the underrwater velocity field, the software uses Airy theory and Stokes weak linear theory for stepper waves with weak non-linearities. The software uses text files to obtain the wave data characteristics as amplitude, period and the fine detail of the grid to calculate the velocities. Th software is written on C without any dependencies to empower its portability on any device/plattaform that can run a GCC compiler.

The software is tought to be able to run on devices as: computers on any plataform (OSX/Linux/Windows), tablets with android, cellphones crossing plattaforms as ARM/X86/PPC as it is a hard only C solution. The software is also made to be faster, future branches will include threading and other processes to improve the code depending on the plattaform it is running on.
