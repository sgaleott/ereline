# totalPipeline.cpp

## MPI processes
- Even processes work on M-arm radiometers, odd processes work on
  S-arm aradiometers.
- This means that any "gather"/"gatherv" operation on gains will mix M
  and S gains together. For instance, having 4 MPI processes will
  split the list of OD to process into *2* chunks, and the result of a
  gainTable::mergeResults() operation will be the following:

                    +-----------+-----------+-----------+-----------+
		MPI process |     #1    |     #2    |     #3    |     #4    |
                    +-----------+-----------+-----------+-----------+
                    |  M gains  |  S gains  |  M gains  |  S gains  |
                    | 1st chunk | 1st chunk | 2nd chunk | 2nd chunk |
                    +-----------+-----------+-----------+-----------+

# ahfInfos.cpp

## ahfInfos::selectProcessRange
- To understand how this function works, read the part about "MPI
  processes" above. The function returns a list containing the number
  of pointing periods that every odd/even chunk must read.
- Here are some outputs (9 ODs):

	    sizeMPI = 2
        detectorIdsSize = 2
        pids_per_od = [ 10 14 9 11 11 8 13 17 12 ]
        nIdsRange (result of ahfInfos::selectProcessRange) = [ 105 ]

        sizeMPI = 4
        detectorIdsSize = 2
        pids_per_od = [ 10 14 9 11 11 8 13 17 12 ]
        nIdsRange (result of ahfInfos::selectProcessRange) = [ 55 50 ]

        sizeMPI = 6
        detectorIdsSize = 2
        pids_per_od = [ 10 14 9 11 11 8 13 17 12 ]
        nIdsRange (result of ahfInfos::selectProcessRange) = [ 44 49 12 ]

        sizeMPI = 8
        detectorIdsSize = 2
        pids_per_od = [ 10 14 9 11 11 8 13 17 12 ]
        nIdsRange (result of ahfInfos::selectProcessRange) = [ 33 30 42 0 ]

        sizeMPI = 10
        detectorIdsSize = 2
        pids_per_od = [ 10 14 9 11 11 8 13 17 12 ]
        nIdsRange (result of ahfInfos::selectProcessRange) = [ 24 20 19 30 12 ]

  Thus, for instance, if we take sizeMPI = 6, we have that each MPI
  job will process a given subset of pointing periods for just *one* arm:

        +-------+--------+--------+--------+
		|       |  4 OD  |  4 OD  |  1 OD  |
		|       | 44 p.p | 49 p.p | 12 p.p | nIdsRange
        +-------+--------+--------+--------+
		| M arm |   #1   |   #3   |   #5   |
		| S arm |   #2   |   #4   |   #6   |
        +-------+--------+--------+--------+

- Note that the implementation of this function does not always
  produce optimal results: the last pair of MPI jobs seems to always
  have much less work than the others. In The case of 9 ODs and 8 MPI
  processes, we have that processes #7 and #8 do not do any work at
  all:

        +-------+--------+--------+--------+--------+
		|       |  3 OD  |  3 OD  |  3 OD  |  0 OD  |
		|       | 33 p.p | 30 p.p | 42 p.p |  0 p.p | nIdsRange
        +-------+--------+--------+--------+--------+
		| M arm |   #1   |   #3   |   #5   |   #7   |
		| S arm |   #2   |   #4   |   #6   |   #8   |
        +-------+--------+--------+--------+--------+

  A better subdivision would have been

        +-------+--------+--------+--------+--------+
		|       |  3 OD  |  2 OD  |  2 OD  |  2 OD  |
		|       | 33 p.p | 22 p.p | 21 p.p | 29 p.p | nIdsRange
        +-------+--------+--------+--------+--------+
		| M arm |   #1   |   #3   |   #5   |   #7   |
		| S arm |   #2   |   #4   |   #6   |   #8   |
        +-------+--------+--------+--------+--------+

- The calls to "floor" and "ceil" are useless, as they are done on
  variables with an "int" types. Removing these calls (and all the
  typecasts) has no impact on the results. (This has been verified.)

## ahfInfos::copyProcessRange
- The baroque use of MPI::COMM_WORLD.Allreduce is to mimick a
  "broadcast" operation. Basically, the value of "nIdsRange" is copied
  from the root process to the others. (Compare this with the main()
  in totalPipeline.cpp, where ahfInfos::selectProcessRange is called
  only in the root process).
- Basically, this method is useless. Just call
  "ahfInfos::selectProcessRange" in every MPI process. (This has been
  verified.)

## ahfInfos::loadAndGetPidDouble
- This function loads a given list.AHF_infos object (thus encompassing
  one OD), it initializes the member variable "firstPointing", and
  finally it returns a list of pointing period IDs as converted to
  floating-pointing numbers (double).
- Some conversion for the "pointID_unique" column is needed, as this
  column stores data as strings. However, It is not clear why such
  numbers are converted to double, as whenever this method is called,
  the elements of the resulting vector are *always* converted to
  integers.
- As a reference, the code that produces the UCDS converts the strings
  in "pointID_unique" into 32-bit integers.

# dipoleFit.cpp
- The members "dipole" and "data" contain the binned dipole signal and the
  binned data (voltage) as a function of the pixel.
- "dipole" and "data" do *not* use the same normalization: the
  "dipole" member contains the *sum* of all the passes over each
  pixel, while "data" has been normalized by the number of hits.

# gainTable.cpp

## gainTable::selectRadiometerGains
- The precondition for this function to work is that gainTable
  contains the M/S gains mixed by gainTable::mergeResults. The outcome
  of the function is to filter the gains so that at the end only the M
  *or* S gains are kept in gainTable (depending on the parameter
  "detectorIdIdx").
- The name "detectorIdIdx" used for the first parameter is misleading:
  this is the *arm* number (0 stands for "M", 1 for "S").
- There are three for cycles in the implementation, but two of them
  can be replaced with calls to std::accumulate and std::insert
  (easier to read). (This has been verified). Also, there are a lot of
  indexes (idx, intIdx, offsetIdx, startPoint) when just two are
  enough.

## gainTable::mergeResults
- The purpose of the function is to collect all the "partial" gain
  tables from each MPI process and merge them so that the object will
  contain the "full" gain table. Using MPI_Gatherv would have
  resulted in much clearer code.
- The code

        int * lengths = new int[sizeMPI];
        fill (lengths,lengths+sizeMPI,0);
        for (int jj=0; jj<sizeMPI; ++jj)
          {
            if (jj==rankMPI)
        	lengths[jj]=length;
            else
        	lengths[jj]=0.0;
          }

  is completely equivalent to

        int * lengths = new int[sizeMPI];
        fill (lengths,lengths+sizeMPI,0);
        lengths[rankMPI] = length;

  which is in turn equivalent to

        std::vector<int> lengths(sizeMPI);
		lengths[rankMPI] = length;

  (the latter also avoid a memory leak).
