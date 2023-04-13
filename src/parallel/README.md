# Welcome to the parallel execution of SpeedPPI

Steps 2 (MSA creation) and 3 (structure prediction) benefit from parallel execution (especially step 3).
\
Using the scripts here, you can evaluate the PPIs in parallel on a computational cluster. If the time you set is too short and the jobs are cancelled - don't worry - the run will continue where it left off automatically.

1. First run step 1 in create_ppi_[some_vs_some/all_vs_all].sh.


2. Modify and run hhblits_parallel.sh.
For the some-vs-some mode this should be run twice - once for each list.

3. Modify and run alphafold_[some_vs_some/all_vs_all]_parallel.sh.

4. Continue with the final steps in create_ppi_[some_vs_some/all_vs_all].sh.
