OLD=v2.0.2
NEW=v2.0.3

download : $(OLD)-job01 $(NEW)-job01

clean :
	rm -rf $(OLD)-* $(NEW)-*

$(OLD)-% :
	rm -rf $@
	mkdir -p $@
	scp -r gamma:profile/$@/log $@/log
	scp -r gamma:profile/$@/mapping_stats.txt $@
	scp -r gamma:profile/map-counts/$(OLD).readnum_to_count $@

$(NEW)-% :
	rm -rf $@
	mkdir -p $@
	scp -r gamma:profile/$@/log $@/log
	scp -r gamma:profile/$@/mapping_stats.txt $@
	scp -r gamma:profile/map-counts/$(NEW).readnum_to_count $@

heatmap_table: v2.0.2-job01/v2.0.2.readnum_to_count v2.0.3-job01/v2.0.3.readnum_to_count
	python mapping_count_dist.py $^

heatmap_table_2: v2.0.2-job01/small v2.0.3-job01/small
	python mapping_count_dist.py $^

countalns : countalns.c
	gcc $^ -o $@
