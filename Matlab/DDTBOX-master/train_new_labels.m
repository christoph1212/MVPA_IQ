function new_labels_train = train_new_labels(labels_train, perm_order)


    for no = 1:size(labels_train, 2)
        new_labels_train(1, no) = labels_train(1, perm_order(no));
    end % of for no
    
end