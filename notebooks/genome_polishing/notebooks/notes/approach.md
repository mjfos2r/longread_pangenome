# Systematic Approach for Manual Plasmid Verification and Genome Polishing

We mainly have agreement between classification methods for most methods, however, there are some particularly interesting discrepancies that are likely due to assembly specific limitations.
Specifically, we need to validate the lack of lp56 calls in several longread assemblies and manually verify the calls across the board.
This is going to be done manually.
I have QUAST output for each assembly. I also have BLAST alignment results for each assembly.
We also need to justify the omission of the strange fusion plasmids from the DB. Jake said that this needs to be checked and justified. My original reason was that we were making assertions regarding correlation/co-presence with confidence that I felt was not justified in the calls.
e.g. if we have a contig that is 2600bp and it hits as a fusion plasmid, that is not enough justification to make further claims. It should be regarded as a fragment and claims made from it's classification should be tempered as such.

Basically me and Ben have been tasked with polishing up and finalizing these genomes. This is going to entail manually validating each contig against the classification that it is being called as. The contigs that have full agreement will not be of any issue, however, we still must visually check them to make sure that things are truly lining up.
The contigs that have discrepant calls will require far more effort to validate.
Fortunately I have already prepared a table of those calls.
The general workflow should be as follows.

This spreadsheet will be made into a google sheet, shared between me and Ben.
The assemblies are accessible between both of us, I can subset those assemblies with the problematic calls.
Or rather, separate out the annotations for those problematic calls. I need to drag over the contig renaming scheme to line up with the plasmid calls like I did for the sql database I threw together last week. I already have that so just need to marry that with the current plasmid call table, rather, just run the same lambda func on the table same as before.

Anyhow, the approach is to be as follows:

1. We have a table of discrepant calls and sequencing method differences. Use this to pull the contigs that have the problematic calls.

2. ....

Actually, first, let's