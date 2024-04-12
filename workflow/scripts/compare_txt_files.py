def count_common_srids(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        srids1 = set(f1.read().split())
        srids2 = set(f2.read().split())
    
    common_srids = srids1.intersection(srids2)
    return len(common_srids)

file1 = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_all_eid.txt'
file2 = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/sample_eid_list.txt'
common_srids_count = count_common_srids(file1, file2)

print(f"Number of common SRIDs in both files: {common_srids_count}")