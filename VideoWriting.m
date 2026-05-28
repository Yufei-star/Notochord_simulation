% Make a video
function []=VideoWriting(L,N,Rt,tdata,R_all,...
    P_all,Pi_all,PE_all,PiE_all,J_loss_all,...
    Pr,tr,FileName,CbOri,FontSetup)

writerObj = VideoWriter(strcat(FileName,'.mp4'),'MPEG-4'); % Name it.
writerObj.FrameRate = 20; % How many frames per second.
open(writerObj);

% Video parameters
n_time_video = 100;
video_time_points = linspace(min(tdata),max(tdata),n_time_video);
temp = abs(video_time_points-tdata');
[~,video_time_idx] = min(temp);

%Input
CLIM=[min(P_all(1,:)), max(P_all(1,:));
    min(Pi_all(1,:)), max(Pi_all(1,:));
    0,4]; % limits for colorbar

for i=1:n_time_video
    disp(['Video_Write_<',FileName,'>_progress:',num2str(i/n_time_video*100),'%'])
    P = P_all(video_time_idx(i),:) * Pr;
    Pi = Pi_all(video_time_idx(i),:) * Pr;
    PE = PE_all(video_time_idx(i),:) * Pr;
    J_loss = J_loss_all(video_time_idx(i),:);
    R = R_all(video_time_idx(i),:);
    plot_notochord_video(L,N,R,Rt,tdata(video_time_idx(i)),P,Pi,PE,[tdata(1:video_time_idx(i))*tr;PiE_all(1:video_time_idx(i),1)'],J_loss/3,tr,CLIM,CbOri,FontSetup);
    drawnow;
    hold off;
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end

disp(['Video successfully written!|frame#:',num2str(n_time_video),'|Frame rate: ',num2str(writerObj.FrameRate)])

end