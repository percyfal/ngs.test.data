ngs.test.data
=============

Test data set generation for ngs


## Installation ## 

Note that you *must* use the *develop* flag.

```bash
python setup.py develop
```

### Issues ###

#### biopython and numpy ####

**biopython** pulls in **numpy** as a dependency, and for some reason,
the **numpy** installation fails. 
  
Solution: installing **numpy** manually first may resolve the issue

```bash
pip install numpy
```
