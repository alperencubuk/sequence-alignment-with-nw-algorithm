## Sequence Alignment with NW Algorithm

### Summary:

Implementation of Needleman–Wunsch algorithm. Compare sequences (in fasta files) and list max scored sequences. Multithreaded with OpenMP.

### How To Run:

- gcc SequenceAlignment.c -o seqAlign -fopenmp
- ./seqAlign

### Notes:

- Increase the LIMIT value (constant) of the program will increase the working time.
- You may need to edit THREAD value (constant) for your processor/operation system.

### Approximately Execution Time:

- 1K Sequences ~ 1 min 30 secs
- 5K Sequences ~ 40 mins 30 secs
- 10K Sequences ~ 2 hours 45 mins
- 15K Sequences ~ 6 hours 15 mins

Tested on **Ubuntu 18.04** with **Intel i7 4710HQ** processor.

---

**Alperen Cubuk**