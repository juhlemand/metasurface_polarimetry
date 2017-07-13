function out = check_beam_power()
    global pwr_meter
    global h_rot_stage
    global h_rot_mount
    global h_rot_mount2
    
    global limit_in_beam
    global limit_out_beam
    
    h_rot_stage.SetAbsMovePos(0, limit_in_beam);
    h_rot_stage.MoveAbsolute(0,1);
    
    tic;
    while and(toc<36, or(ismoving(h_rot_mount)==1, ismoving(h_rot_mount2)==1))
       pause(1) 
    end
    
    fprintf(pwr_meter, 'MEAS:POW?');
    out = str2num(fscanf(pwr_meter));
    
    h_rot_stage.SetAbsMovePos(0, limit_out_beam);
    h_rot_stage.MoveAbsolute(0,1);
    

end
