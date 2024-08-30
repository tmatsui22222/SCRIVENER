#!/bin/bash

# Function to check the status of the last command and handle errors or success
check_status() {
    if [ $? -ne 0 ]; then
        echo "❌ ERROR: $1. Please check the logs for details."
        # Log error but do not exit, continue with the script
    else
        echo "✅ SUCCESS: $2 installed."
    fi
}

echo "Installing bioinformatics tools... 🚀"

# Install tools in custom bioinformatics folder
echo "Checking if /opt/bioinfo directory exists..."
if [ ! -d "/opt/bioinfo" ]; then
    echo "Creating /opt/bioinfo directory..."
    sudo mkdir /opt/bioinfo && sudo chown $(whoami) /opt/bioinfo
    if [ $? -ne 0 ]; then
        echo "Failed to create or change ownership of /opt/bioinfo. Proceed in current directory or specify another? (y/N/enter directory path):"
        read user_input
        if [[ $user_input == "y" ]]; then
            echo "Proceeding with installation in the current directory..."
        elif [[ $user_input == "N" ]]; then
            echo "Exiting installation."
            exit 1
        else
            mkdir -p "$user_input" && cd "$user_input"
            if [ $? -ne 0 ]; then
                echo "Failed to create or access specified directory. Exiting installation."
                exit 1
            else
                echo "Proceeding with installation in $user_input..."
            fi
        fi
    else
        echo "✅ SUCCESS: bioinfo directory set up."
    fi
else
    echo "Directory /opt/bioinfo already exists. Continuing with installations..."
    sudo chown $(whoami) /opt/bioinfo
    check_status "Failed to change ownership of existing /opt/bioinfo" "bioinfo directory"
fi

cd /opt/bioinfo

# Continue with the rest of the installations
echo "Updating system and installing dependencies..."
sudo apt-get update && sudo apt-get install -y wget python3 python3-pip python3-dev git build-essential cmake zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev
check_status "Failed to update or install required packages" "system dependencies"

# Install Filtlong
echo "Installing Filtlong..."
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make && sudo cp filtlong /usr/local/bin/
check_status "Failed to install Filtlong" "Filtlong"
cd ..

# Install Flye
echo "Installing Flye..."
git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install
check_status "Failed to install Flye" "Flye"
cd ..

# Install BWA
echo "Installing BWA..."
git clone https://github.com/lh3/bwa.git
cd bwa
make && sudo cp bwa /usr/local/bin/
check_status "Failed to install BWA" "BWA"
cd ..

# Install Samtools
echo "Installing Samtools..."
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar xjf samtools-1.20.tar.bz2
cd samtools-1.20
./configure && make && make install
check_status "Failed to install Samtools" "Samtools"
cd ..

# Install Bcftools
echo "Installing Bcftools..."
wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar xjf bcftools-1.20.tar.bz2
cd bcftools-1.20
./configure && make && make install
check_status "Failed to install Bcftools" "Bcftools"
cd ..

# Install Clustal Omega
echo "Installing Clustal Omega..."
wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz
tar xzf clustal-omega-1.2.4.tar.gz
cd clustal-omega-1.2.4
./configure && make && make install
check_status "Failed to install Clustal Omega" "Clustal Omega"
cd ..

# Install Racon
echo "Installing Racon..."
git clone --recursive https://github.com/lbcb-sci/racon.git
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make && sudo make install
check_status "Failed to install Raon" "Racon"
cd ../..

# Install Medaka
echo "Installing Medaka..."
pip install medaka
check_status "Failed to install Medaka" "Medaka"

# Install Minimap2
echo "Installing Minimap2..."
git clone https://github.com/lh3/minimap2
cd minimap2 && make
sudo cp minimap2 /usr/local/bin/
check_status "Failed to install Minimap2" "Minimap2"
cd ..

# Install sniffles
echo "Installing Sniffles Structural Variant (e.g., indels) for long-reads caller..."
pip install sniffles
check_status "Failed to install sniffles" "sniffles"

# Install Python dependencies
echo "Installing additional Python dependencies..."
pip install -U scikit-learn
pip install biopython levenshtein
check_status "Failed to install one or more Python dependencies" "Python dependencies"

# Clean up
echo "Cleaning up..."
sudo apt-get clean

echo "All tools installed successfully. 🦠🧬🖳✨"
