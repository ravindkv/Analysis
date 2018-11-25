mass=$1
channel=$2
category=$3

echo "INPUT FILE: " $mass
cp limitComputaionCondor_13TeV.py ../
sed -i "s:MASS_POINT_:$mass:g" limitComputaionCondor_13TeV.py
sed -i "s:CHANNEL_:$channel:g" limitComputaionCondor_13TeV.py
sed -i "s:CAT_:$category:g" limitComputaionCondor_13TeV.py

python limitComputaionCondor_13TeV.py
mv ../limitComputaionCondor_13TeV.py .
