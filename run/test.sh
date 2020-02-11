for((i=0;i<$1;i++))
do
    fn=` printf "model_output/run%04d/I211_J211.dat" $i `
    if [ ! -f "${fn}" ]; then echo $i "The file doesn't exist"; fi
    fn=` printf "model_output/run%04d/I2212_J2212.dat" $i `
    if [ ! -f "${fn}" ]; then echo $i "The file doesn't exist"; fi
    fn=` printf "model_output/run%04d/I321_J2212.dat" $i `
    if [ ! -f "${fn}" ]; then echo $i "The file doesn't exist"; fi
    fn=` printf "model_output/run%04d/I321_J321.dat" $i `
    if [ ! -f "${fn}" ]; then echo $i "The file doesn't exist"; fi
done
