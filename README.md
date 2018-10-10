# exSeek
Also see my [website page](https://www.cmwonderland.com/blog/2018/10/10/101-exrna-project/)

The aim of this project is to combine **small cell free RNA sequencing techniques** with **machine learning** to identify promising candidate biomarker for cancers. We establish the pipeline from mapping to feature selection and propose a **robust feature selection method** to identify potential biomarkers.


### mapping 


We use sequential mapping strategy for small RNA data (HCC exRNA and GSE71008)
We use fast mapping order test method to map in thousands of different orders in a very short time and decide the order.
![Markdown](http://i2.tiimg.com/640680/a7331057c61ed2a8.png)

### construct expression matrix

Considering the characteristic of small RNA and exRNA. We think it is better to use domain features instead of full length RNA. So we also do peak calling to construct both full length and peak expression matrix.

### Preprocessing
All the datasets need to be normalized and remove batch effect.

I test several normalization method:

- CPM
- remove piRNA and miRNA
- Use several reference miRNA
- SCNorm
- TMM
- RLE

I also use batch removal method including RUVs and Combat.

###  Robust feature selection method
Our aim is to establish a robust feature selection machine learning model to select very few feature as candidate cancer biomarkers.

- Normalize domain coverage 
- Scale each feature (log) independently (using z-scores, min-max, **robust normalization**)
- Run a classifier (**Random Forest, Logistic Regression, Linear SVM**) to select features based on feature importance. Optimize hyper-parameters by **3-fold cross-validation**.
- Optionally, use a recursive feature elimination (**RFE**) to eliminate features.
- Repeat feature selection for 100 times. Randomly remove 10% samples in each run (**Jackknife resampling**) or do **shuflle split**.
- Select features that are recurrently selected across resampling runs (>50%)
- Refit the classifier on selected features

For detailed illustration of our pipeline please see [**Work Summary**](#work)



## Related
As teaching assistant of [Bioinformatics Basics](https://legacy.gitbook.com/book/lulab/teaching/details). I use part of our work to assign the [final project](https://lulab.gitbooks.io/teaching/content/quiz/quiz_exrna/quiz_exrna_tutorial.html) of this semester's course.  
