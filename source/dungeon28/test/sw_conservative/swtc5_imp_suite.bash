#!/bin/bash
#PBS -A cli106ms
#PBS -q debug
#PBS -N homme_cpu
#PBS -l gres=atlas1
#PBS -l nodes=3
#PBS -j oe -l walltime=1:00:00

# to run, first build executables 

#these are needed to use the GPU using OpenACC and they change frequently
#so check info about the target machine
export CRAY_CUDA_MPS=1
export CRAYPE_LINK_TYPE="dynamic"
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

#these are all similar to swtc5_bd.nl file in repo
export SCRIPT1=swtc5_exp.nl
export SCRIPT2=swtc5_rk.nl
export SCRIPT3=swtc5_imp.nl
export SCRIPT4=swtc5_imp2.nl
export SCRIPT5=swtc5_bd.nl

# executable for implicit code
#GNU version
export EXE=/ccs/home/4ue/trunk/np4/src/swim/swim
#PGI version
export EXE2=/ccs/home/4ue/trunk/pgi_cpu_read/src/swim/swim
#GPU OpenACC version
export EXE3=/ccs/home/4ue/trunk/pgi_gpu_read/src/swim/swim
#run directory
export BASE=/lustre/atlas1/cli106/proj-shared/4ue/test/swtc5_table1

cd $PBS_O_WORKDIR

echo leapfrog GNU
cp -f $SCRIPT1 input.nl
time aprun -n48 $EXE < $BASE/input.nl  &> lf_120s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_lf_120s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_lf_120s

echo RK GNU
cp -f $SCRIPT2 input.nl
time aprun -n48 $EXE < $BASE/input.nl  &> rk_30s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_rk_30s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_rk_30s

echo RK GNU PGI
cp -f $SCRIPT2 input.nl
time aprun -n48 $EXE2 < $BASE/input.nl  &> rk_pgi_30s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_rk_pgi_30s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_rk_pgi_30s

echo Crank-Nicolson GNU
cp -f $SCRIPT3 input.nl
time aprun -n48 $EXE < $BASE/input.nl  &> cn_120s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_cn_120s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_cn_120s

echo Crank-Nicolson BIG TS
cp -f $SCRIPT4 input.nl
time aprun -n48 $EXE < $BASE/input.nl  &> cn_1800s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_cn_1800s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_cn_1800s

echo BDF2 BIG TS
cp -f $SCRIPT5 input.nl
time aprun -n48 $EXE < $BASE/input.nl  &> bd_1800s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_bd_1800s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_bd_1800s

echo Crank-Nicolson BIG TS PGI
cp -f $SCRIPT4 input.nl
time aprun -n48 $EXE2< $BASE/input.nl  &> cn_pgi_cpu_1800s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_cn_pgi_1800s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_cn_pgi_1800s

echo Crank-Nicolson BIG TS PGI GPU
cp -f $SCRIPT4 input.nl
time aprun -n48 $EXE3< $BASE/input.nl  &> cn_pgi_gpu_1800s.out
mv -f HommeSWTime_stats HommeSWTime_stats_ne30np4_cn_pgi_gpu_1800s
mv -f HommeSWTime.00 HommeSWTime.ne30pn4_cn_pgi_gpu_1800s
