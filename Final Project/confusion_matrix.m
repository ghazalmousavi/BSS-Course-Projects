function cm = confusion_matrix(cm_old, true_labels, predicted_labels)
    cm = cm_old;

    for i = 1:length(true_labels)
        cm(true_labels(i), predicted_labels(i)) = cm(true_labels(i), predicted_labels(i)) + 1;
    end
    
end
