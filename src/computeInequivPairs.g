#
# Usage:
#   gap -q src/computeInequivPairs.g
#   gap -q -c 'n:=13; chunks:=6; runChunk:=2; Read("src/computeInequivPairs.g");'
# Args:
#   n: group size; chunks: number of chunks to split work; runChunk: 1-based chunk index.
#

if not IsBound(n) then
    n := 10;
fi;
if not IsBound(chunks) then
    chunks := 1;
fi;
if not IsBound(runChunk) then
    runChunk := 1;
fi;

start_time := NanosecondsSinceEpoch();

G := SymmetricGroup(n);;
CC := ConjugacyClasses(G);
pairs := [];
for i in [1..Length(CC)] do
    for j in [i..Length(CC)] do
        Add(pairs,[CC[i],CC[j]]);
    od;
od;

SizeOfChunks := Int(Floor(Float(Length(pairs)/chunks)));
ChunkIndicies := [1,SizeOfChunks+1..SizeOfChunks*chunks+1];
if runChunk = chunks then
    pairs := pairs{[ChunkIndicies[runChunk]..Length(pairs)]};
else
    pairs := pairs{[ChunkIndicies[runChunk]..ChunkIndicies[runChunk+1]-1]};
fi;
uniquePairs := [];
for pair in pairs do
    rep1 := ListPerm(Representative(pair[1]), n);
    rep2 := Representative(pair[2]);
    c1 := Centralizer(pair[1]);
    c2 := Centralizer(pair[2]);
    doubleCosetRepsAndSizes := DoubleCosetRepsAndSizes(G, c1, c2);
    for repAndSize in doubleCosetRepsAndSizes do
        rep := repAndSize[1];
        newPair := ListPerm(rep2^rep, n);
        Add(uniquePairs,[rep1,newPair]);
    od;
od;

path := Directory("pairs_data");
# Output file name would "Pairs-11-6-2" means DS11 Pairings, with this being the second chunk of 6
filename := Concatenation("Pairs-",String(n),"-",String(chunks),"-",String(runChunk),".txt");
file := Filename(path, filename);
output := OutputTextFile(file , false);;

for pair in uniquePairs do
    AppendTo(output, pair, "\n");
od;

CloseStream(output);
Display(Concatenation(["Saved to: ", file]));
s := Concatenation(["n = ", String(n), ". Pairs = ", String(Length(uniquePairs))]);
Display(s);
end_time := Float((NanosecondsSinceEpoch() - start_time)/1*10^(-9));
Display(Concatenation(["Time taken (min): ", String(Float(end_time/60))]));

quit;
