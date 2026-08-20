# simpleNMRjeolTools
interface between JEOL Jason NMR files and simpleNMR

## Running the program

Run the program from the directory where the file ***simpleNMRjeolTools.py*** is from the comand line

```python simpleNMRjeolTools.py```



## Running the program for the first time

The server side program checks if the user is registered each time it runs and so since this will be the first time you have run it you have to give an email address. We don't save much on the server side, just an email address and the id of the computer and count the number of times the user runs the program.

Once the user has registered and accepted the terms the user has to run the program again and everything should work. I haven't checked this from this program as I am already registered! Hopefully it will work!

## Notes

- I have only tested the program on two examples and found some problems, I think i have sorted them out, but you never know.
- One of the problems I found was the creation of the molecule by extracting the information from the h5 file. It worked for the first molecule, but not the second, now it works for both, but who knows if it works for every molecule.
- The other problem I found when running the program on the second example was due to an indexing problem, hopefully I have fixed that for other examples too.
