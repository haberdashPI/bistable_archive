
#!/bin/sh

result_dir=$1
group=${2:-3}
name=${3:-"reinterpret"}

proj_dir="projects/bistable/src"

# S=101 # start from where N=100 leaves off
# N=10
# N=2000 # start with just 10 jobs for now, and see how that goes.

cd ${proj_dir}
GIT_HASH=`git rev-parse HEAD`
cd

files=(${result_dir}/data/*.jld2)
mkdir -p ${result_dir}/logs/$name

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq 0 $group ${#files[@]}`; do
  file_group=${files[@]:i:group}
  echo "sbatch ${proj_dir}/run_reinterpret.sh \
	  $file_group --git_hash $GIT_HASH \
	  --params ${result_dir}/params.jld2 \
	  --settings ${proj_dir}/settings.toml \
	  --folder $name \
	  -l ${result_dir}/logs/$name/result_${i}.log"
done

