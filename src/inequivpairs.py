def make_inequiv_pairs(n):
    """ Convert from gap permutation pairs output to txt format
        gap format is assumed to be stored in Pairs-n-1-1.txt
    """
    path = "data10/"
    fname = f"Pairs-{n}-1-1.txt"
    with open(f"{path}{fname}", "r") as f:
        lst = f.readlines()
    inequiv_parser_pairs(lst, f"{path}New_{fname}", n)

def inequiv_parser_pairs(lst, out_file, n, save_chunk=100000):
    """ parse gap txt output in lst
        make new format
    """
    pairs_lst = []
    count = 0
    with open(out_file, "w") as f:
        i = 0
        while i < len(lst):
            r = lst[i]
            row_count = 0
            # cycles of first perm in pair
            middle = r[1:-1]
            left_paren = middle.find("[")
            right_paren = middle.find("],")
            first_perm = middle[left_paren:right_paren+1]
            if first_perm.find(".") < 0:
                first_perm = first_perm.split()
                first_lst = [int(i.strip(",")) for i in first_perm if i.strip(",").isdigit()]
            else:
                first_lst = [i+1 for i in range(n)] # is identity

            # cycles of second perm in pair
            sec_half = middle[right_paren+1:]
            sec_left_paren = sec_half.find("[")
            if sec_left_paren < 0: # on a half line
                i += 1
                r = lst[i]
                middle = r[1:-1]
                sec_left_paren = middle.find("[")
                sec_half = middle
            sec_right_paren = sec_half.rfind("]")
            sec_perm = sec_half[sec_left_paren:sec_right_paren+1]
            if sec_perm.find(".") < 0:
                sec_perm = sec_perm.split()
                sec_lst = [int(i.strip(",")) for i in sec_perm if i.strip(",").isdigit()]
            else:
                sec_lst = [i+1 for i in range(n)] # is identity

            pairs_lst.append((first_lst, sec_lst))
            count += 1
            if count >= save_chunk:
                # save memory by writing in chunks
                f.writelines((str(pair)[1:-1]+"\n" for pair in pairs_lst))
                del pairs_lst
                pairs_lst = []
                count = 0
            i+=1

        f.writelines((str(pair)[1:-1]+"\n" for pair in pairs_lst))
    del pairs_lst
    print("done saving")

def read_inequiv_pairs(n):
    """ read New_Pairs-n-1-1.txt and return pairs of permutations
    """
    with open(f"data10/New_Pairs-{n}-1-1.txt", "r") as f:
        lst = f.readlines()
    perms = []
    for r in lst:
        middle = r[0:-1]
        left_paren = middle.find("[")
        right_paren = middle.find("],")
        first_perm = middle[left_paren+1:right_paren]
        first_perm = first_perm.replace(",", " ").split()
        first_lst = [int(i) for i in first_perm]

        sec_half = middle[right_paren+1:]
        sec_left_paren = sec_half.find("[")
        sec_right_paren = sec_half.rfind("]")
        sec_perm = sec_half[sec_left_paren+1:sec_right_paren]
        sec_perm = sec_perm.replace(",", " ").split()
        sec_lst = [int(i) for i in sec_perm]
        perms.append((first_lst, sec_lst))
    return perms

