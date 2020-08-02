# UGER Basics

#### Accessing UGER:

`ssh kgarg@login.broadinstitute.org`

**Interactive Sessions:**

Request memory on a per core basis \(total = mem \* \# cores\)

`ish -l h_vmem=8G -binding linear:4 -pe smp 4` 

\(4 cores and 32 GB\) 

_defaults_: 1GB mem; 2 hr runtime; 1 core \([link](https://broad.service-now.com/sp?id=kb_article_view&sys_kb_id=63354be4137f0b00449fb86f3244b049)\)

#### Storgage Space:

`/broad/hptmp/`

**Interactive Sessions:**

\*\*\*\*

