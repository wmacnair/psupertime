Hi Dario!

For the moment, I haven't built the package, so you need to use devtools to run it. The documentation should all be there though.

Do your git clone / pull, then to load the functions you need to do: 
```R
devtools::load_all('/path/to/psupertime')
```

Let me know if this causes any problems.

Cheers
Will