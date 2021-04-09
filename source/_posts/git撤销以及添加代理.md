---
title: git撤销操作，添加代理以及wget设置代理
date: 2021-04-09 11:14:38
tags: 转载
categories: linux
description: 这些东西都能搜到，但也容易忘记，记录在此方便查阅。
---

## git撤销操作

### 没有add

```bash
git status
git checkout . # 撤销修改。
git checkout <somefile> # 加文件名撤销该文件的修改。
```

### 没有commit

```bash
git status
git reset HEAD # 撤销上一次add.
git reset HEAD <somefile> # HEAD后加文件名，撤销该文件。
```

### 没有push

```bash
git log # 得到你需要回退一次提交的commit id
git reset --hard <commit_id>  # 回到其中你想要的某个版本
git reset --hard HEAD^  # 回到最新的一次提交
git reset HEAD^  # 此时代码保留，回到 git add 之前
```

### 已经push，希望回滚到上一个版本

1. 通过git reset直接删除指定的commit

```bash
git log # 得到你需要回退一次提交的commit id
git reset --hard <commit_id>
git push origin HEAD --force # 强制提交一次，之前提交就从远程仓库删除
```

2. 通过git revert是用一次新的commit来回滚之前的commit

```bash
git log # 得到你需要回退一次提交的commit id
git revert <commit_id>  # 撤销指定的版本，撤销也会作为一次提交进行保存
```

git revert是用一次新的commit来回滚之前的commit，此次提交之前的commit都会被保留；
git reset是回到某次提交，提交及之前的commit都会被保留，但是此commit id之后的修改都会被删除。

## git添加代理

```bash
git config --global --add remote.origin.proxy <代理地址> # 添加代理
git config --global --add remote.origin.proxy 127.0.0.1:8889
git config --global --unset remote.origin.proxy # 取消代理
```

## wget设置代理

添加变量即可

```bash
export HTTP_PROXY=127.0.0.1:8889
export HTTPS_PROXY=127.0.0.1:8889
# 取消
unset HTTP_PROXY
unset HTTPS_PROXY
```