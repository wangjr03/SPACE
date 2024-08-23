Refer to DLPFC_data_analysis.md for step-by-step analysis code based on the step 1 and step 2 results from the SPACE algorithm. 

For larger datasets, the first two steps of SPACE are typically performed using a high-performance computing cluster by breaking the data into smaller gene chunks. This makes the code runs faster and smoother. Here is the description of the process: 

Step 1: run code code_Step1.R using the .sb file step1.sb.
Based on how small the gene chunks(10 genes in the example code) are and the total no of genes in the dataset, the number N_ARRAY is decided. N_ARRAY=ceiling(total no of genes in dataset/chunk size). In this example of sample 6, N_ARRAY=ceiling(4865/10)=487. 

combining results from the chunk output:

```{r}
sample_i=6
data=matrix(0,nrow=1,ncol=12)
colnames(data)=c("X","GSP1","COS1","GSP2","COS2","GSP3","COS3","GSP4","COS4","GSP5","COS5","combined")
for(iter in 1:N_ARRAY){
    r=read.csv(paste0("dataset/out_step1_sample",sample_i,"_",iter,".csv"))
    data=rbind(data,r)
}

data=data[-1,]
write.csv(data,paste0("dataset/Results/Sample_",sample_i,"_result_step1.csv"))

p_adj3=p.adjust(data[,12], method = "BY")
a3=which(p_adj3<0.05)
length(a3)
```

Step 2: run code code_Step2.R using the same .sb file step1.sb. Just make sure to change the name of the Rscript and N_ARRAY. For this step, N_array=ceiling(a3/chunk size) 
Note: Step 2 uses step 1 output as input. If you prefer to use other methods to detect SVG in step 1. Simply omit the step 1 and do step 2 with slight adjustment in code. 

combining results from the output of step 2:
```{r}
data=matrix(0,nrow=1,ncol=12)
colnames(data)=c("X","GSP1","COS1","GSP2","COS2","GSP3","COS3","GSP4","COS4","GSP5","COS5","combined")
for(i in 1:N_ARRAY){
    r=read.csv(paste0("dataset/Results/out_step2_sample_",sample_i,"_",i,".csv"))
    data=rbind(data,r)
}
data=data[-1,]

write.csv(data,paste0("dataset/Results/sample_",sample_i,"_result_step2.csv"))


data=matrix(0,nrow=1,ncol=2)
colnames(data)=c("X","V1")
for(i in 1:N_ARRAY){
    r=read.csv(paste0("dataset/Results/list_step2_sample_",sample_i,"_",i,".csv"))
    data=rbind(data,r)
}
list1=data[-1,]

list1[,1]=1:dim(list1)[1] #adjustments for computing chunk wise
write.csv(list1,paste0("dataset/Results/sample_",sample_i,"_result_step2_list1.csv"))
```
Now these outputs can be used in subsequent steps provided in DLPFC_data_analysis.md

