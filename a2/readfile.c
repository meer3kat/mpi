int read_file(char *name, int** pp){
	//return the number of number read in the file and pointer *pp will point to the first
    FILE* f;
    /* if need to only enter file name for data in input_data/ folder
    char dest[100] = "input_data/";
    strcat(dest, name);
    printf("%s", dest);
    f = fopen(dest,"rb");*/
    f = fopen(name, "r");
    int count = 0;
    if(f){
        //printf("file opened! \n");
        fseek(f, 0, SEEK_END);
        fseek(f, 0, SEEK_SET);
        int *p = NULL;
        int tmp;
        int n;
        fscanf(f,"%d ",&n);
        p = (int*)malloc(n*sizeof(int));
        for(int i=0;i<n;i++){
        	fscanf(f,"%d ", &p[i]);
        }
        // while(fscanf(f,"%d ", &tmp)>0){
        // 	p[count] = tmp;
        // 	count++;
        // 	if(count==initial_size){
        // 		initial_size = initial_size*2;
        // 		p = (int*)realloc(p,2*initial_size*sizeof(int));
        // 	}
        // }
        // p = (int*)realloc(p,count*sizeof(int));
        *pp = &p[0];
        fclose(f);
        return n;
    }
    else
    	return 0;
}