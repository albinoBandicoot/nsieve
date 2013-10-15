.file "matrow_ops_32.s"
.text
.globl asm_xor_row
.globl bitscan 

asm_xor_row:
	movl 	12(%esp), %edx
	movl	8(%esp), %ecx
	movl	4(%esp), %eax
LOOP:
	cmpl	$1, %edx
	jle	DONE
	movdqu	(%eax), %xmm0
	movdqu	(%ecx), %xmm1
	pxor	%xmm0, %xmm1
	movdqu	%xmm1, (%eax)
	addl	$16, %eax
	addl	$16, %ecx
	subl	$2, %edx
	jmp	LOOP
DONE:
	cmpl	$1, %edx
	jne	RETURN
	movl	(%ecx), %edx
	xorl	%edx, (%eax)
	movl	4(%ecx), %edx
	xorl	%edx, 4(%eax)
	ret
RETURN:
	ret
	
bitscan:
	cmpl	$0, 4(%esp)
	je	HIGH_ORDER
	bsr	4(%esp), %eax
	ret
	
HIGH_ORDER:
	bsr	8(%esp), %eax
	addl	$32, %eax
	ret
	
