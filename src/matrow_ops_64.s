.file "matrow_ops_64.s"
.text
.globl _xor_row
.globl _bitscan 

_xor_row:
LOOP:
	cmpl	$1, %edx
	jle	DONE
	movdqu	(%rdi), %xmm0
	movdqu	(%rsi), %xmm1
	pxor	%xmm0, %xmm1
	movdqu	%xmm1, (%rdi)
	add	$16, %rdi
	add	$16, %rsi
	decl	%edx
	decl	%edx
	jmp	LOOP
DONE:
	cmpl	$1, %edx
	jne	RETURN
	mov	(%rsi), %rdx
	xor	%rdx, (%rdi)
RETURN:
	ret
	
_bitscan:
	bsf	%edi, %eax
	ret
	
