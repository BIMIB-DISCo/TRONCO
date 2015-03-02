#### hypothesis.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Add a new hypothesis by creating a new causal event and adding it to the dateset
"hypothesis.add" <-
  function( data, label.formula, lifted.formula,  ... ) {

    label.effect = list(...);
    #print(label.effect)
    
    # save the needed data structures
    if(!is.null(data$genotypes) && !is.null(data$annotations)) {
      dataset = data$genotypes;
      annotations = data$annotations;
    }
    else {
      dataset = NULL;
      annotations = NULL;
    }
    if(!is.null(data$hypotheses)) {
      hypotheses = data$hypotheses;
    }
    else {
      hypotheses = NA;
    }
    
    # add the hypothesis only if all the inputs are correctly provided
    if(!is.null(dataset) && !is.null(annotations)) {
      # the Boolean functions look for a global variable named lifting.dataset
      # if there are already global variables named as the ones used here, make the backup of them
      do.roll.back.lifting.dataset = FALSE;
      do.roll.back.lifting.annotations = FALSE;
      do.roll.back.lifting.edges = FALSE;
      
      # I need a global variable to save the dataset of the lifted formula
      # if there is already a global variable named lifting.dataset, make the backup of it
      if(exists("lifting.dataset")) {
        roll.back.lifting.dataset = lifting.dataset;
        do.roll.back.lifting.dataset = TRUE;
      }
      assign("lifting.dataset",dataset,envir=.GlobalEnv);
      
      # I need a global variable to save the annotations of the lifted formula
      # if there is already a global variable named lifting.annotations, make the backup of it
      if(exists("lifting.annotations")) {
        roll.back.lifting.annotations = lifting.annotations;
        do.roll.back.lifting.annotations = TRUE;
      }
      assign("lifting.annotations",annotations,envir=.GlobalEnv);
      
      # I need a global variable to save the edges of the lifted formula
      # if there is already a global variable named lifting.edges, make the backup of it
      
      do.roll.back.lifting.dataset = FALSE; # <- ???? why????
      
      if(exists("lifting.edges")) {
        roll.back.lifting.edges = lifting.edges;
        do.roll.back.lifting.edges = TRUE;
      }
      assign("lifting.edges",NULL,envir=.GlobalEnv);
      
      #print('lifted.formula')
      #print(lifted.formula)
      
      ## test
      #lifted.formula = eval(lifted.formula)
      
      
      # save the lifted dataset and its hypotheses for the current formula
      curr_formula = lifted.formula$formula;
      curr_hypotheses = lifted.formula$hypotheses;
      


      #print('cur formula')
      #print(curr_formula)
      #print('cur_hypo')
      #print(curr_hypotheses)
      
      # save the edges of the lifted formula
      hstructure = lifting.edges;
      
      # roll back to the previous value of the global variable lifting.dataset if any or remove it
      if(do.roll.back.lifting.dataset) {
        assign("lifting.dataset",roll.back.lifting.dataset,envir=.GlobalEnv);
      }
      else {
        rm(lifting.dataset,pos=".GlobalEnv");
      }
      
      # roll back to the previous value of the global variable lifting.annotations if any or remove it
      if(do.roll.back.lifting.annotations) {
        assign("lifting.annotations",roll.back.lifting.annotations,envir=.GlobalEnv);
      }
      else {
        rm(lifting.annotations,pos=".GlobalEnv");
      }
      # roll back to the previous value of the global variable lifting.edges if any or remove it
      if(do.roll.back.lifting.edges) {
        assign("lifting.edges",roll.back.lifting.edges,envir=.GlobalEnv);
      }
      else {
        rm(lifting.edges,pos=".GlobalEnv");
      }
      # set the hypotheses number
      if(!is.na(hypotheses[1])) {
        num.hypotheses = hypotheses$num.hypotheses;
      }
      else {
        num.hypotheses = 0;
      }
      # * is a special label.effect which indicates to use all the events as effects for this formula
      is.to.all.effects = FALSE;
            
      if(label.effect[[1]][1]=="*") {
        label.effect = colnames(dataset)[1:(length(colnames(dataset))-num.hypotheses)];
        #any event can not be both causes and effects for the formula to be well-formed
        label.effect = list(label.effect[-which((label.effect%in%unlist(curr_hypotheses$llist)))]);
        is.to.all.effects = TRUE;
        
        if(length(label.effect)==0) 
        {
          stop(paste("[ERR] Missing list of effects to test or wildcard \'*\'.", sep=''));
        }
      }
      # check the formula to be well-formed
      all.col.nums = vector();
      if(length(label.effect)==0) {
        stop(paste("[ERR] Missing list of effects or wildcard \'*\'.", sep=''));
      }
      else {
        #check the effects of the formula to be well-formed
        for (i in 1:length(label.effect)) {
          curr.label.effect = label.effect[[i]];
          if(is.to.all.effects==FALSE) {
            col.num = -1;
            if(length(curr.label.effect)==1) {
              event.map = emap(c(curr.label.effect,"*"),dataset,annotations);
              col.num = event.map$col.num;
              events.name = event.map$events.name;
            }
            else if(length(curr.label.effect)==2) {
              event.map = emap(curr.label.effect,dataset,annotations);
              col.num = event.map$col.num;
              events.name = event.map$events.name;
            }
          }
          else {
            col.num = which(colnames(dataset)%in%curr.label.effect);
            if(length(col.num)==0) {
              col.num = -1;
            }
            events.name = curr.label.effect;
          }
          #check the effect to be a valid event
          if(col.num[1]==-1) {            
            stop(paste("[ERR] Unknown gene among effects: \"", curr.label.effect,
                       "\".",sep=''));
          }
          all.col.nums = append(all.col.nums,col.num);
          #check the formula to be well-formed
          #if the effect is in the formula, the formula is not well-formed
          if(length(which(unlist(curr_hypotheses$llist)%in%events.name))>0) {
                  stop(paste("[ERR] Bad forme formula, event \"", curr.label.effect,
                       "\" yields a loop.",,sep=''));          
            }
        }
      }
      
      # look for duplicated effects in the formula
      if(anyDuplicated(all.col.nums)>0) 
        {
        stop(paste("[ERR] Bad formed formula, duplicated events ", 
                   paste(label.effect[duplicated(label.effect)], collapse=', ', sep=''),
                   "within effects.", sep=''));          
        }
      #check that the we are not duplicating any name by adding the new hypothesis
      if(length(which(colnames(dataset)==label.formula))>0) 
      {
        stop(paste("[ERR] This hypothesis already exists.", sep=''));
      }
      #add the hypothesis to the dataset
      dataset = cbind(dataset,curr_formula);		
      #check that the formula is valid according to Suppes' theory
      #structure to compute the observed and observed joint probabilities
      pair.count <- array(0, dim=c(ncol(dataset),ncol(dataset)));
      #compute the probabilities on the dataset
      for(i in 1:ncol(dataset)) {
        for(j in 1:ncol(dataset)) {
          val1 = dataset[ ,i];
          val2 = dataset[ ,j];
          pair.count[i,j] = (t(val1) %*% val2);
        }
      }
      #marginal.probs is an array of the observed marginal probabilities
      marginal.probs <- array(as.matrix(diag(pair.count)/nrow(dataset)),dim=c(ncol(dataset),1));
      #joint.probs is an array of the observed joint probabilities
      joint.probs <- as.matrix(pair.count/nrow(dataset));
      #check that the probability of the formula is in (0,1)
      if(marginal.probs[ncol(dataset)]==0 || marginal.probs[ncol(dataset)]==1) 
        {        
        stop(paste("[ERR] Formula has marginal probability ", marginal.probs[ncol(dataset)], 
                   ", but should be in (0,1).", sep=''));
      }
      #check that the formula does not duplicate any existing column
      i = ncol(dataset);
      for(j in 1:ncol(dataset)) {
        #if the edge is valid, i.e., not self cause
        if(i!=j) {
          #if the two considered events are not distinguishable
          if((joint.probs[i,j]/marginal.probs[i])==1 && (joint.probs[i,j]/marginal.probs[j])==1) 
          {
            stop(paste("[ERR] The formula duplicates event (", paste(as.events(data)[j, ], collapse=', ', sep=''), 
                       ").", sep=''));
          }
        }
      }
      #now I can finally add the hypothesis
      colnames(dataset)[ncol(dataset)] = label.formula;
      if(is.na(hypotheses[1])) {
        hypotheses = list();
      }
      hypotheses$num.hypotheses = num.hypotheses + 1;
      #create the list of added hypotheses
      if(length(hypotheses$hlist)==0) {
        hypotheses$hlist = vector();
      }
      #add the new hypothesis to the list
      for (i in 1:length(label.effect)) {
        curr.label.effect = label.effect[[i]];
        if(is.to.all.effects==FALSE) {
          if(length(curr.label.effect)==1) {
            event.map = emap(c(curr.label.effect,"*"),dataset,annotations);
            col.num = event.map$col.num;
          }
          else if(length(curr.label.effect)==2) {
            event.map = emap(curr.label.effect,dataset,annotations);
            col.num = event.map$col.num;
          }
        }
        else {
          col.num = which(colnames(dataset)%in%curr.label.effect);
          if(length(col.num)==0) {
            col.num = -1;
          }
        }
        for (j in 1:length(col.num)) {
          hypotheses$hlist = rbind(hypotheses$hlist,t(c(colnames(dataset)[ncol(dataset)],colnames(dataset)[col.num[j]])));
        }
        if(is.null(colnames(hypotheses$hlist))) {
          colnames(hypotheses$hlist) = c("cause","effect");
        }
      }
      #create the list of hypotheses' structures
      if(length(hypotheses$hstructure)==0) {
        hypotheses$hstructure = new.env(hash=TRUE,parent=emptyenv());
      }
      hypotheses$hstructure[[label.formula]] = get.lifted.formula(hstructure);
      #add the new hypothesis in the annotations
      annotations = rbind(data$annotations,c("Hypothesis", label.formula));
      rownames(annotations)[nrow(annotations)] = label.formula;
      #add the color of the type "Hypothesis" is not already defined
      if(any(rownames(data$types)=="Hypothesis")==FALSE) {
        types = rbind(data$types, 'slateblue');
        rownames(types)[nrow(types)] = "Hypothesis";
        data$types = types;
      }
      #return the new data as result
      data$genotypes = dataset;
      data$hypotheses = hypotheses;
      data$annotations = annotations;
      return(data);
    }
    else {
      stop("[ERR] Missing dataset or formula.");
    }
    return(NA);
  }

#### end of file -- hypothesis.add.R
