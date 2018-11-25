echo "\section{Cat-2: Integral of inclusive Mjj from b-jet pT bins(10 tables)}"
./runMe.sh 90 MuChannel bJetPtCat

echo "\newpage"
echo "\section{Cat-4: Integral of Mjj from exclusive charm-categories (3 tables)}"
./runMe.sh 90 MuChannel CTagCat

echo "\newpage"
echo "\section{Cat-5: Integral of Mjj from exclusive charm-categories in bins of b-jet pT (30 tables)}"
./runMe.sh 90 MuChannel BothCat

rm HplusShapes_*
rm datacard_csbar*
