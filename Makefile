test:
	python sv_multimatching.py -i utils/test_datasets/variant_test.bed.gz -r utils/test_datasets/ref_test_1.bed.gz
	python sv_multimatching.py -i utils/test_datasets/variant_test.bed.gz -r utils/test_datasets/ref_test_2.bed.gz
	python sv_multimatching.py -i utils/test_datasets/variant_test.bed.gz -r utils/test_datasets/ref_test_3.bed.gz
	python sv_multimatching.py -i utils/test_datasets/variant_test.bed.gz -r utils/test_datasets/ref_test_4.bed.gz
	python sv_multimatching.py -i utils/test_datasets/variant_test.bed.gz -r utils/test_datasets/ref_test_5.bed.gz
	python sv_multimatching.py -i utils/test_datasets/variant_test.bed.gz -r utils/test_datasets/ref_test_6.bed.gz
