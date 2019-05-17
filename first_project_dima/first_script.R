matt <- GSE122020_series_matrix

write.csv(matt, file='gse122020_matrix_IDasInd.csv',row.names = 1)
read.csv('gse122020_matrix_IDasInd.csv')

write.csv(matt, file='gse122020_matrix_IDasCol.csv')
read.csv('gse122020_matrix_IDasCol.csv')


