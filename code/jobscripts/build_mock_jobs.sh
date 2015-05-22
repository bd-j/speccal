ncpu=256
nwalker=`expr $ncpu - 1`
nwalker=`expr $nwalker \* 2` 
niter=1024
datadir=$PROJECTS/speccal/data/ggclib/mocks/miles/

files=$(ls ${datadir}ggc_mock.u0*.pkl | xargs -n 1 basename)

params=$(ls ${datadir}ggc_mock.u0*pkl | xargs -n 1 basename | sed -e "s/ggc_mock.\(.*\)\.pkl/\1/")
for p in $params
do
    #tmp = ${f#*.}
    #p = ${tmp
    sed -e "s/\${nwalker}/$nwalker/" -e "s/\${niter}/$niter/" \
	-e "s/\${params}/$p/" -e "s/\${ncpu}/$ncpu/" \
	stampede_template_mocks.sh > stampede_mocks_$p.sh
    #echo $p
done
