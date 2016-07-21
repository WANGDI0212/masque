#!/bin/sh

name=AlienTrimmer
javac $(pwd)/$name.java ;
if [ ! -e $(pwd)/$name.class ]; then exit; fi;

mkdir $(pwd)/JARBUILDER/
cp $(pwd)/$name.class $(pwd)/JARBUILDER/ ;

mkdir $(pwd)/JARBUILDER/META-INF/ ;
echo "Main-Class: $name" > $(pwd)/JARBUILDER/META-INF/MANIFEST.MF ;

cd ./JARBUILDER/ ;
jar cvfm $name.jar ./META-INF/MANIFEST.MF . ;
cd ../ ;

cp JARBUILDER/$name.jar . ;

rm -r $(pwd)/JARBUILDER/ ; 
rm $(pwd)/$name.class ;
