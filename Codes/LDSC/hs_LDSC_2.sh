sumstats_files=(
    "2hGlu_2021"
    "EPI_ALL_2018"
    "LDLC_UKB"
    "AD_2019"
    "EPI_DR_2014"
    "MDD_2018"
    "ADHD_2019"
    "EPI_FD_2018"
    "NEUR_2018"
    "AF_UKB"
    "EPI_GEN_2018"
    "PD_2019"
    "ALS_2018"
    "EPI_ILAE_2014"
    "RTC_UKB"
    "AN_2017"
    "FG_2021"
    "SCZ_2021"
    "ANX_CC_2016"
    "FI_2021"
    "SD_2019"
    "ANX_FS_2016"
    "FTD_2017"
    "SWB_2016"
    "ASD_2017"
    "Glu_UKB"
    "systolic_bp"
    "BD_2019"
    "HbA1c_2021"
    "systolic"
    "BMI_UKB"
    "HDLC_UKB"
    "T2D"
    "CAD"
    "HV_2017"
    "UKB_460K_body_HEIGHTz"
    "CHD_UKB"
    "IHD"
    "UKB_460K_cov_EDU_YEARS"
    "diastolic"
    "IHD_UKB"
    "UKB_460K_pigment_HAIR"
    "DS_2016"
    "INS_2019"
    "UKB_460K_pigment_SKIN"
    "EA_2018"
    "INT_2018"
    "UKBB_COVID19"
)

TXT=.txt
SUMSTAT=.sumstats.gz
for i in *.txt;do
    for sumstats_file in "${sumstats_files[@]}"; do
#basefile_name=$(basename $i .txt)
basefile_name=$basename$i$TXT
#sumstats_name=$(basename $sumstats_file .sumstats.gz)
sumstats_name=$basename$sumstats_file$SUMSTAT

./ldsc.py \
--h2 /rds/general/user/zw4419/home/LDSC/SumStats/$sumstats_file \
--w-ld-chr /rds/general/user/zw4419/home/LDSC/REF/weights_hm3_no_hla/weights. \
--ref-ld-chr /rds/general/user/zw4419/home/LDSC/Annot_ALL/$basefile_name.,/rds/general/user/zw4419/home/LDSC/REF/1000G_EUR_Phase3_baselineNEW/baseline. \
--overlap-annot \
--frqfile-chr /rds/general/user/zw4419/home/LDSC/REF/1000G_Phase3_frq/1000G.EUR.QC. \
--out $basefile_name.$sumstats_name \
--print-coefficients
mv $basefile_name.$sumstats_name.results /rds/general/user/zw4419/home/LDSC/LDSC_Results_FATSL
    done
done