;;; export.el --- Export orgmode docs to RST for Sphinx -*- lexical-binding: t -*-

(require 'org)
(require 'ox-rst)

(defun radsex-export-all ()
  "Export all orgmode files under docs/orgmode/ to RST under docs/source/."
  (let ((org-dir (expand-file-name "orgmode" default-directory))
        (rst-dir (expand-file-name "source" default-directory)))
    (dolist (org-file (directory-files-recursively org-dir "\\.org$"))
      (let* ((relative (file-relative-name org-file org-dir))
             (rst-file (expand-file-name
                        (concat (file-name-sans-extension relative) ".rst")
                        rst-dir))
             (rst-parent (file-name-directory rst-file)))
        (unless (file-directory-p rst-parent)
          (make-directory rst-parent t))
        (with-current-buffer (find-file-noselect org-file)
          (org-export-to-file 'rst rst-file nil nil nil nil)
          (kill-buffer))
        (message "Exported %s -> %s" org-file rst-file)))))

(radsex-export-all)
