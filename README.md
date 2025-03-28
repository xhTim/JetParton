This is the framework to couple the jet parton shower and parton cascade

# Install LHAPDF
1. **Download the LHAPDF Package**  
   Download the LHAPDF compressed file
2. mkdir /Tool/LHAPDF_Lib
3. cd /path/to/extracted-LHAPDF-directory
4. ./configure --prefix=/Tool/LHAPDF_Lib --disable-python
5. make -j8

# Install PDF Files
1. wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF31_nnlo_as_0118.tar.gz
2. tar -xzvf NNPDF31_nnlo_as_0018.tar.gz
3. mv NNPDF31_nnlo_as_0018 /Tool/LHAPDF_Lib/share/LHAPDF/

# Compile pythia8 with LHAPDF
1. download pythia8.310
2. After extracting Pythia, enter the resulting directory.
3. ./configure --prefix=/Tool/pythia8310_install --with-lhapdf6=/Tool/LHAPDF_install
4. make -j8
5. make install
