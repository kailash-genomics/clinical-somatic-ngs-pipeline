#!/usr/bin/env bash
set -e

echo "[INFO] Starting pipeline setup"

BASE_DIR="/content/drive/MyDrive/Somatic_Demo_Analysis"
TOOLS_DIR="/content/tools"

mkdir -p "$TOOLS_DIR"
cd "$TOOLS_DIR"

#####################################
# SYSTEM & JAVA
#####################################
apt-get update -qq
apt-get install -y openjdk-17-jre wget unzip

#####################################
# Core bioinformatics tools
#####################################
command -v fastqc   &>/dev/null || apt-get install -y fastqc
command -v fastp    &>/dev/null || apt-get install -y fastp
command -v bwa      &>/dev/null || apt-get install -y bwa
command -v samtools &>/dev/null || apt-get install -y samtools
command -v bcftools &>/dev/null || apt-get install -y bcftools

#####################################
# htslib (bgzip + tabix)  ← THIS IS NEW
#####################################
command -v bgzip &>/dev/null || apt-get install -y tabix

#####################################
# GATK 4
#####################################
if ! command -v gatk &>/dev/null; then
  wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
  unzip -q gatk-4.5.0.0.zip
  ln -s "$TOOLS_DIR/gatk-4.5.0.0/gatk" /usr/local/bin/gatk
fi

#####################################
# SnpEff (JAR already in Drive)
#####################################
if ! command -v snpEff &>/dev/null; then
  ln -s /content/drive/MyDrive/snpEff/snpEff.jar /usr/local/bin/snpEff.jar
  cat << 'WRAP' > /usr/local/bin/snpEff
#!/usr/bin/env bash
java -Xmx4g -jar /usr/local/bin/snpEff.jar "$@"
WRAP
  chmod +x /usr/local/bin/snpEff
fi

#####################################
# Required directories
#####################################
mkdir -p \
  "$BASE_DIR/alignment" \
  "$BASE_DIR/metrics" \
  "$BASE_DIR/logs" \
  "$BASE_DIR/reports"

echo "[INFO] Setup completed successfully"
