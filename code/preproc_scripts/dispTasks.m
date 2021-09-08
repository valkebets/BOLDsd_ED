function dispTasks(tasksTodo, tasksDone)
% display tasks for preprocessing script
    global CLASSICSEGMENT;
    global NEWSEGMENT;

    disp('The following tasks have already been done:')
    if tasksDone.realign, disp('-Realignment');  end
    if isfield(tasksDone, 'QC') && tasksDone.QC, disp('-Quality Control'); end
    if tasksDone.coregister, disp('-Coregistration'); end
    if tasksDone.coregReslice, disp('-coregReslice'); end
    if tasksDone.smooth, disp('-smooth'); end
    if tasksDone.segment == NEWSEGMENT, disp('-New Segmentation'); end
    if tasksDone.segment == CLASSICSEGMENT, disp('-Classic Segmentation'); end
    if tasksDone.label, disp('-Labeling'); end
    
    disp('The following tasks have to be done:')
    if tasksTodo.realign, disp('-Realignment');  end
    if tasksTodo.QC, disp('-Quality Control'); end
    if tasksTodo.coregister, disp('-Coregistration'); end
    if tasksTodo.coregReslice, disp('-coregReslice'); end
    if tasksTodo.smooth, disp('-smooth'); end
    if tasksTodo.segment == NEWSEGMENT, disp('-New Segmentation'); end
    if tasksTodo.segment == CLASSICSEGMENT, disp('-Classic Segmentation'); end
    if tasksTodo.label, disp('-Labeling'); end
return
    