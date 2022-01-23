# rose tutorial

ROSE(Rank Ordering of Super-Enhancers) is a tool for detecting super-enhancers.

```
# @ host
$ docker pull hattyoriiiiiii/rose:1.0.0
$ docker run -it --rm -v $PWD:/home/hattori/work -w /home/hattori hattyoriiiiiii/rose:1.0.0

# @ container
$ cd rose
$ python2 ROSE_main.py -g HG18 -i ./data/HG18_MM1S_MED1.gff -r ./data/MM1S_MED1.hg18.bwt.sorted.bam -c ./data/MM1S_WCE.hg18.bwt.sorted.bam -o example_test -s 12500 -t 2500
```