import os
#reorder unaligned xy sections
def reorder_sections(exp_folder,starts,stops,excludes):
    ol = []
    numseg = len(starts)
    print numseg
    print excludes
    for i in range(0,numseg):
        if starts[i]>stops[i]:
            ol = ol + range(starts[i],stops[i]-1,-1)
        else:
            ol = ol + range(starts[i],stops[i]+1)
    utfolders = os.listdir(exp_folder + "unaligned/")
    #print utfolders
    print len(ol)-len(excludes)
    for utfolder in utfolders:
        if os.path.isdir(exp_folder + "unaligned/" + utfolder):
	    ims = sorted(os.listdir(exp_folder + "unaligned/" + utfolder))
            print exp_folder + "unaligned/" + utfolder
            reord = 0
            for j in ol:
            	oldfile = (exp_folder + "unaligned/" + utfolder + "/" + ims[j-1])
                #print oldfile
                if not j in set(excludes):
                    reord = reord + 1
                    newfile = (exp_folder + "unaligned/" +
                               utfolder + "/" + "%04d" % reord + ".tif")
                    #print newfile
                    os.rename(oldfile, newfile)
                else:
                    os.remove(oldfile)
    return
