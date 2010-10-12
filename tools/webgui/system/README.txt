Compilation of the wrappers.

1. Compile static version of the library - find main Makefile
 located in the trunk directory and type:
 
 make static
 
 After this step you new file (libamtrack.a or libamtrack.lib) should
 appear in lib directory.
 
 2. Compile all wrappers - go to each directory containing wrapper
 code and type also:
 
 make static
 
 New .exe file should be created.
 
 3. Remember to copy all .exe files to the wrapper directory on the server. 