# equaids 1.2.1 -- the copy submitted to the SSC archive

The same code as equaids 1.2.0 on GitHub; only the version stamps differ. It is
here to test the SSC package before the archive publishes it. As on SSC, the
package file lists every file with an `f` line: `net install` installs the
programs, the help and the dialog only; the example data and the other
ancillary files are copied by `net get` to the current folder, never to the
system directories.

```stata
net install equaids, from("https://raw.githubusercontent.com/aabbdd12/equaids/main/ssc/1.2.1") replace
net get equaids, from("https://raw.githubusercontent.com/aabbdd12/equaids/main/ssc/1.2.1") replace
```

The examples of the help read the data from the current folder, else from the
SSC archive, else from GitHub. To go back to the GitHub version:

```stata
net install equaids, from("https://raw.githubusercontent.com/aabbdd12/equaids/v1.2.0") replace
```
