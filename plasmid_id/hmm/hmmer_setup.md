# Let's type some pfam32s and classify some plasmids!
ok so I have a multifasta testor set from Ira Schwarz containing 36 plasmid pfam32 types.
I have used muscle to generate an MSA file in the .afa format
## Let's use some muscle to align all of the testors!

```{bash}
muscle -align Minimal-PF32-testor-set-March2021.faa -output pf32_ref_v1.afa
```

## Let's install hmmer + easel
We need hmmer to actually build and use the HMM, but before building, we need to convert our `.afa` to `.sto`.

I have installed hmmer with brew, but the brew installation did not come with easel.
Thus I shall compile it from source using the following installation instructions from the hmmer documentation: (after uninstalling brew's hmmer...)
### Installation/compillation instructions for hmmer
 - to download and build the current source code release:

```{bash}
   % wget http://eddylab.org/software/hmmer/hmmer.tar.gz
   % tar zxf hmmer.tar.gz
   % cd hmmer-3.4
   % ./configure --prefix /your/install/path   # replace /your/install/path with what you want, obv
   % make
   % make check                                # optional: run automated tests
   % make install                              # optional: install HMMER programs, man pages
   % (cd easel; make install)                  # optional: install Easel tools
```

I am installing it to /usr/local (*THIS REQUIRES ROOT*)

```{bash}
./configure --prefix /usr/local
make
make check
sudo make install
cd easel
sudo make install
```

everything is now in path, including the easel applets!

Now to convert my `.afa` to `.sto`!

```{bash}
esl-reformat stockholm pf32_ref_v1.afa > pf32_ref_v1.sto
```

## Let's build our HMM!
okay so now that we have installed hmmer+easel, we have performed MSA, converted our MSA to stockholm, let's build our HMM.

continuing from the documentation, the commands required to build an HMM profile are:

```{bash}
hmmbuild pf32_ref_v1.hmm pf32_ref_v1.sto
```

and now we can take a look at our shiny new hmm!
>[! HMM profile]-
>```{hmm}
>HMMER3/f [3.4 | Aug 2023]
>NAME  pf32_ref_v1
>LENG  249
>ALPH  amino
>RF    no
>MM    no
>CONS  yes
>CS    no
>MAP   yes
>DATE  Tue Mar 12 14:07:20 2024
>NSEQ  36
>EFFN  1.270020
>CKSUM 1085671701
>STATS LOCAL MSV      -10.9539  0.70293
>STATS LOCAL VITERBI  -11.5681  0.70293
>STATS LOCAL FORWARD   -5.4179  0.70293
>HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
>            m->m     m->i     m->d     i->m     i->i     d->m     d->d
>  COMPO   2.67773  4.52230  2.96278  2.67789  3.26494  3.18427  3.83155  2.51392  2.54134  2.31911  3.57177  2.97561  3.49926  3.12563  3.15728  2.62455  2.85485  2.61460  4.94747  3.36536
>          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
>          0.01981  4.32714  5.04949  0.61958  0.77255  0.00000        *
>      1   3.40170  4.88936  4.58335  4.23641  3.50427  4.20324  4.86643  2.73093  3.97515  2.04889  0.71219  4.45028  4.68635  4.33661  4.12165  3.78681  3.72677  2.78631  5.39003  4.18399      1 m - - -
>          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
>          0.07677  4.32714  2.80199  0.61958  0.77255  0.48576  0.95510
>      2   2.84498  5.60241  0.91629  2.24979  4.94342  3.34556  3.93699  4.47677  2.98099  3.98799  4.85123  2.69842  3.98893  2.91190  3.58041  2.97354  3.38835  4.04013  6.14808  4.63325      2 d - - -
>          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
>          0.05064  4.27133  3.34062  0.61958  0.77255  0.56096  0.84551
>          ............... (many more lines)
>```

okay, so now lets run our multifasta file containing all of the parA sequences extracted from my ONT assemblies (see pfam_typing.ipynb)

but first, let's run the multifasta we built the HMM from through to see if everything remains groovy.

```{bash}
hmmsearch pf32_ref_v1.hmm Minimal-PF32-testor-set-March2021.faa
```

I think it's groovy? Perhaps I need to dupe the file but replace the headers with numbers to make sure that what I'm seeing in the output is what I think it is.

[ update ] okay that is finished. Each record has been assigned the header of "TEST_#" and the corresponding renaming headers has been written to a table for lookup.

ok now let's rerun that hmmsearch but this time on our test file.

[*ALSO I HAVE SORTED THE DIRECTORY. ALL SEQS ARE IN plasmid_id/seqs, all HMMS are in /plasmid_id/hmm*]

```{bash}
hmmsearch hmm/pf32_ref_v1.hmm seqs/pf32_ref_test_set.faa > hmm_out/pf32_ref_v1_test.out
```

>[!pf32_ref_v1_test.out partial output]-
># hmmsearch :: search profile(s) against a sequence database
># HMMER 3.4 (Aug 2023); http://hmmer.org/
># Copyright (C) 2023 Howard Hughes Medical Institute.
># Freely distributed under the BSD open source license.
># - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
># query HMM file:                  hmm/pf32_ref_v1.hmm
># target sequence database:        seqs/pf32_ref_v1_test_set.faa
># - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
>
>Query:       pf32_ref_v1  [M=249]
>Scores for complete sequences (score includes all domains):
>   --- full sequence ---   --- best 1 domain ---    -#dom-
>    E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
>    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
>   4.2e-111  361.2  18.2   5.1e-111  361.0  18.2    1.0  1  TEST_35
>   2.7e-110  358.6  15.9     3e-110  358.4  15.9    1.0  1  TEST_31
>   4.3e-109  354.6  17.4   5.4e-109  354.3  17.4    1.0  1  TEST_33
>   9.1e-108  350.3  17.4     1e-107  350.2  17.4    1.0  1  TEST_26
>     8e-107  347.2  22.0   9.4e-107  347.0  22.0    1.0  1  TEST_30
>     5e-106  344.6  14.3     6e-106  344.4  14.3    1.0  1  TEST_28
>   1.3e-105  343.2  18.2   1.6e-105  342.9  18.2    1.0  1  TEST_21
>   1.4e-104  339.9  14.5   1.5e-104  339.7  14.5    1.0  1  TEST_25
>     1e-103  337.0  13.5   1.2e-103  336.9  13.5    1.0  1  TEST_29
>   3.7e-103  335.2  18.9   4.1e-103  335.1  18.9    1.0  1  TEST_23
>   1.1e-102  333.7  18.8   1.2e-102  333.6  18.8    1.0  1  TEST_13
>     7e-102  331.0  21.9   7.8e-102  330.9  21.9    1.0  1  TEST_10
>   2.3e-101  329.3  13.9   2.9e-101  329.0  13.9    1.0  1  TEST_18
>     3e-101  329.0  20.6   3.3e-101  328.8  20.6    1.0  1  TEST_19
>   5.1e-101  328.2  19.8   5.6e-101  328.1  19.8    1.0  1  TEST_36
>   1.6e-100  326.6  19.9   1.8e-100  326.4  19.9    1.0  1  TEST_5
>   1.8e-100  326.4  15.7     2e-100  326.2  15.7    1.0  1  TEST_1
>   7.2e-100  324.4  17.1     8e-100  324.3  17.1    1.0  1  TEST_7
>    1.6e-99  323.3  22.6    1.8e-99  323.2  22.6    1.0  1  TEST_27
>    2.6e-99  322.6  15.4    3.1e-99  322.4  15.4    1.0  1  TEST_34
>      3e-99  322.4  21.0    3.3e-99  322.3  21.0    1.0  1  TEST_22
>    5.1e-99  321.7  15.3    6.2e-99  321.4  15.3    1.0  1  TEST_17
>    4.6e-98  318.5  17.2    5.1e-98  318.4  17.2    1.0  1  TEST_15
>    5.2e-98  318.3  18.3    5.8e-98  318.2  18.3    1.0  1  TEST_20
>    6.7e-98  318.0  13.7    8.2e-98  317.7  13.7    1.0  1  TEST_8
>    1.1e-97  317.3  18.9    1.2e-97  317.1  18.9    1.0  1  TEST_9
>      2e-97  316.4  19.4    2.2e-97  316.3  19.4    1.0  1  TEST_11
>    1.1e-96  314.0  22.4    1.2e-96  313.9  22.4    1.0  1  TEST_4
>    3.7e-96  312.3  14.2    4.1e-96  312.2  14.2    1.0  1  TEST_12
>    1.8e-95  310.0  19.8      2e-95  309.9  19.8    1.0  1  TEST_6
>    2.2e-94  306.5  15.7    2.7e-94  306.2  15.7    1.0  1  TEST_3
>      9e-94  304.5  15.5      1e-93  304.3  15.5    1.0  1  TEST_32
>    1.2e-93  304.1  19.1    1.3e-93  303.9  19.1    1.0  1  TEST_24
>      6e-91  295.2  14.1    6.7e-91  295.1  14.1    1.0  1  TEST_14
>    4.8e-90  292.3  10.3    5.3e-90  292.1  10.3    1.0  1  TEST_2
>    3.3e-59  191.2   1.4      4e-59  191.0   1.4    1.1  1  TEST_16
>
>
>Domain annotation for each sequence (and alignments):
>>> TEST_35
>   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
> ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
>   1 !  361.0  18.2  5.1e-111  5.1e-111       1     247 [.       1     245 [.       1     249 [. 0.99
>
>  Alignments for each domain:
>  == domain 1  score: 361.0 bits;  conditional E-value: 5.1e-111
>  pf32_ref_v1   1 mdrkkpkiitiasikGGvGKstlsiifatlLskkkKvLLiDlDpQasltSYflkkikek.vdieklniyevlkekldinkviikindnldlipsylsl 97
>                  md+kkpkiitiasikGGvGKst++iifa+lL++k+KvLLiD+D+Qas+tSYf+k+i ++ ++i ++niy+vlkekldin++ii+i+dnldlipsylsl
>      TEST_35   1 MDNKKPKIITIASIKGGVGKSTSAIIFANLLAQKYKVLLIDIDTQASTTSYFYKEIASQkINIVSKNIYRVLKEKLDINNAIINIKDNLDLIPSYLSL 98
>                  9********************************************************999************************************** PP
>
>  pf32_ref_v1  98 ekfnkesieykelllkkklellkenydYiiiDtsPsldlllknaLvisdyviiPvtaekwsvEslelleeaieklekkkikifiletqfkknnthkel 195
>                  +kf++e i++kel+lk++l +lk++ydYiiiDt+Psld++l+naL++s+++i+P+taekw+vEslelle++i+kl++ ki+ifil+t+fkknnthkel
>      TEST_35  99 HKFSSEFIPLKELRLKDNLIFLKQDYDYIIIDTNPSLDFTLSNALMTSNCIIVPMTAEKWAVESLELLEFYIKKLRI-KIPIFILITRFKKNNTHKEL 195
>                  *****************************************************************************.******************** PP
>
>  pf32_ref_v1 196 lkllkkkyeekflgkiheredlkkliaeneefdlkedYykeykevLekilkk 247
>                  lk++++k  + flg++heredl+++i+ n+efd+++dY++eyke+++k++++
>      TEST_35 196 LKHVRSK--KGFLGIVHEREDLNRKITGNDEFDMTKDYIEEYKEAVSKFFDM 245
>                  *******..***************************************9975 PP
>
>>> TEST_31
>   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
> ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
>   1 !  358.4  15.9    3e-110    3e-110       1     246 [.       1     244 [.       1     246 [] 0.99
>
>  Alignments for each domain:
>  == domain 1  score: 358.4 bits;  conditional E-value: 3e-110
>  pf32_ref_v1   1 mdrkkpkiitiasikGGvGKstlsiifatlLskkkKvLLiDlDpQasltSYflkkikek.vdieklniyevlkekldinkviikindnldlipsylsl 97
>                  md+kkpkiitiasikGGvGKst+siifatlL++k+KvLLiDlD+Qas+tSYf+kki+++ +d+ ++niy+vlk++ld+n++i++i++nldlipsy++l
>      TEST_31   1 MDTKKPKIITIASIKGGVGKSTSSIIFATLLAQKYKVLLIDLDTQASTTSYFCKKIENQkIDLVNKNIYRVLKDTLDVNNAIVNIKENLDLIPSYITL 98
>                  9************************************************************************************************* PP
>
>  pf32_ref_v1  98 ekfnkesieykelllkkklellkenydYiiiDtsPsldlllknaLvisdyviiPvtaekwsvEslelleeaieklekkkikifiletqfkknnthkel 195
>                  +kf++e i+++el+lk++l +lk++ydYi++Dt+Psld++l+naL++s++vi+P+taekw++Esl+lle++ie+l++ ki+if+l+t+fkknnthkel
>      TEST_31  99 HKFSNEFIPHQELRLKDSLIFLKQDYDYIVVDTNPSLDFTLSNALITSNCVIVPMTAEKWAIESLDLLEFHIENLKI-KIPIFLLVTRFKKNNTHKEL 195
>                  *****************************************************************************.******************** PP
>
>  pf32_ref_v1 196 lkllkkkyeekflgkiheredlkkliaeneefdlkedYykeykevLekilk 246
>                  lk+++++  e+flg iheredl+k+ia n+ef++++dY++eyke+L+k+++
>      TEST_31 196 LKYVESR--ERFLGFIHEREDLNKKIAGNNEFNMDKDYINEYKEALSKFFE 244
>                  *******..****************************************96 PP
>
>>> TEST_33  ......
> ....

ok so I do not know how to go from that to the specific alignment within the hmm. Let's step away from the bike shed and keep reading the documentation.

## Okay, let's make a *profile database* instead
The current method is working as intended, however the above approach is more for identifying domains. My HMM above is little more than a limited pf32 set. If it aligns to that HMM, it's a pf32. Yes, I already know that. I want it to differentiate between the different pf32 subtypes.
This means that I need to make an *HMM profile database* where each plasmid type is it's own HMM, and the entire database is used for classification?

this means I need to do the following:

1. download all known plasmids from NCBI
2. using the annotations for each plasmid, separate out the parA genes and then group by source plasmid. (make multifasta for each plasmid subtype)
3. MSA on each multifasta
4. convert afas to stockholms
5. hmmbuild on each stockholm
6. cat all hmms together
7. hmmpress the build hmm for hmmscan
8. repeat hmmscan of test.

### Building the big HMM db

I currently have all of the plasmid *fastas* downloaded, I need the genbank files. After downloading all of those, I'll iterate through each one, and then separate out the parA genes and put them in a multifasta file for each plasmid type.

~~~```{bash}
efetch -something something- > something something
```~~~

okay, so actually I'll just make a notebook for this. Refer to: "multi_hmm_construction.ipynb"
I will be using biopython + efetch to pull the genbanks for each and save them to their own folder based on which plasmid it is.
