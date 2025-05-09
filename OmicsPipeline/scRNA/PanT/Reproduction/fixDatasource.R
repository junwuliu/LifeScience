metainfo = mca@meta.data
for (i in 1:dim(metainfo)[1]){
	list_a = unlist(strsplit(metainfo[i,'datasource'],split="_"))
	if (list_a[1] %in% c("GSE156728","GSM4743237","GSM4743231","GSM4743199")){
		metainfo[i,"newdatasource"] = paste(list_a[2],list_a[1],list_a[3],sep="_")
		metainfo[i,"CancerType"] = list_a[2]
	}else{
		metainfo[i,"newdatasource"] = metainfo[i,"datasource"]
		metainfo[i,"CancerType"] = list_a[1]
	}
}
