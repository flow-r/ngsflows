# Contributing

Each file in the R folder represents a tool. In it we have a R function preferably of the same name, or a similar easy name.

These are some of the practices we follow in-house. We feel using these makes stitching custom pipelines using a set of modules quite easy. Consider this a check-list of a few ideas and a work in progress.

## A note on module functions


1. should accept minimum of **two inputs**, 
    - a input file etc, depends on the module. Flexible
    - samplename (is used to append a column to the flowmat)
  ```
  x
  samplename = get_opts("samplename")
  ```

2. should always return a list arguments:
    - **flowmat** (required)   : contains all the commands to run
    - **outfiles** (recommended): could be used as an input to other tools
  ```
  return(list(outfiles = mergedbam, flowmat = flowmat))
  ```

3. can define all other default arguments such as paths to tools etc. in a seperate conf (tab-delimited) file.
  - Then use `get_opts("param")` to use their value.

 ```
 ## Example conf file:
 cat my.conf
 bwa_exe	/apps/bwa/bin/bwa
 ```

4. should use `check_args()` to make sure none of the default parameters are null. 

 ```{r}
 ## check_args(), checks ALL the arguments of the function, 
 ## and throws a error. use ?check_args for more details.
 get_opts("my_new_tool")
 ```

## Example

```{r picard_merge, echo=TRUE, comment=""}
picard_merge <- function(x,
        samplename = get_opts("samplename"),
        mergedbam,
        java_exe = get_opts("java_exe"),
        java_mem = get_opts("java_mem"),
        java_tmp = get_opts("java_tmp"),
        picard_jar = get_opts("picard_jar")){
	
  ## Make sure all args have a value (not null)
  ## If a variable was not defined in a conf. file get_opts, will return NULL
  check_args()  
  
  bam_list = paste("INPUT=", x, sep = "", collapse = " ")
  ## create a named list of commands
  cmds = list(merge = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s MergeSamFiles %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
  java_exe, java_mem, java_tmp, picard_jar, bam_list, mergedbam))
  
  ## Create a flowmat
  flowmat = to_flowmat(cmds, samplename)
  
  ## return a list, flowmat AND outfiles
  return(list(outfiles = mergedbam, flowmat = flowmat))
}
```


<!-- Here are a few things to note regarding naming a function and what it should do:

- it is *recommended* to have a lower case function name, seperated by `_`. And in general we try to follow Advanced R's [style guide](http://adv-r.had.co.nz/Style.html).
- Each function has:
  - a few input files,
  - a few paths (to files and tools)
  - and default parameters
- Each function return a list, with elements:
  - outfiles: a list/character vector of output file names
  - flowmat: a data.frame, with a few extra attributes -->
  
