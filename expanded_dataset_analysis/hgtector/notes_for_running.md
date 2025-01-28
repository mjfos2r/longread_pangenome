# Running HGTector2

Running an analysis for Amine!
He has given me the file: `LD83_GCs_min1.faa` which I am to run against the standard database.

The only issue is that I need to reinstall the database :^)

Fortunately I have it stored in the bucket at `gs://jel-seq-databases/db/hgtector`

It's also like, 150gb.

Anyway let me clear up some room for this on the VM so I can actually run it.
> Compressing the old AVA results and scuttling to the bucket.

anyway, this will probably be a tomorrow morning problem if I'm being honest.
Anyway we're gonna give it a good try anyway.

Groovy, downloading everything from the bucket. Will compress, checksum, and reupload as a tarball.


## Running Amine's Data
I'm gonna run Amine's data through this tool after I validate the database.

Database validation should be easily done with the following command:
Per the documentation:
>Break and resume
>Should any of the download steps be interrupted by e.g., a network failure, one can resume the downloading process by
>re-executing the same command. **The program will skip the already downloaded files in this new run**
>In some instances, one may need to manually remove the last file from the failed run (because that file may be corrupt), before re-running the program.
>**If one wants to overwrite downloaded files (e.g., upgrading), add --overwrite to the command.**

Ok groovy. let's see what we've got here.

***

Command:
```bash
hgtector database -c bacteria -s 1 -r genus\
                --reference --compile diamond\
                --threads 30\
                -o /home/mf019/db/hgtector/
```
That did indeed work.

***

Now let's run this file for amine!

Command:
```bash
hgtector search -i LD83_GCs_min1.faa -o search_dir \
-m diamond -p 30 \
-d /home/mf019/db/hgtector/diamond/db.dmnd \
-t /home/mf019/db/hgtector/taxdump
```

and now to analyze:

Command:
```bash
hgtector analyze -i search_dir \
-o analyze_dir -t /home/mf019/db/hgtector/taxdump \
--donor-name
```
>>Reran the above with the addition of `--donor-name`
***

Cool! Finished running. Now to compress and send to Amine.

Command:
```bash
tar -cvf - hgtector/ | pigz -p 30 >amine_hgtector_20250122.tar.gz \
    && md5sum amine_hgtector_20250122.tar.gz > amine_hgtector_20250122.tar.gz.md5 \
    && md5sum -c amine_hgtector_20250122.tar.gz.md5
```

#>>{MJF - 2025-Jan-23}<<#
