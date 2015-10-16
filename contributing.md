Contributing
------------------------------------

Contributions to this package is very appreciated.

Each file in the R folder represents a tool. In it we have a R function prefereably of the same name, or a similar easy name.

Here are a few things to note regarding naming a function and what it should do:

- it is *recommended* to have a lower case function name, seperated by `_`. And in general we try to follow Advanced R's [style guide](http://adv-r.had.co.nz/Style.html).
- Each function has:
  - a few input files,
  - a few paths (to files and tools)
  - and default parameters
- Each function return a list, with elements:
  - outfiles: a list/character vector of output file names
  - flowmat: a data.frame, with a few extra attributes

  
  